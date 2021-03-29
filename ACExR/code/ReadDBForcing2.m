function ReadDBForcing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function = ReadDBForcing;
%
% Reads boundary forcing data for 'loch' from 
% database sources e.g. daily wind speed, daily rainfall/
%    river flow, coastal density profile
% Assumes daily values of data over 1 year.
%
% Usage:    LochData contains the catalogue data.
%           Bdata is a structured array containing 
%                   boundary data.
%
% Phil Gillibrand
% 14/3/2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('Reading boundary forcing data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set simulation length
Ndays = Param.Ndays + Param.ndays_spin_up;
ndays_spin_up = Param.ndays_spin_up;
Param.Ndays = Ndays;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read freshwater input data.
% if Bdata.Qf_switch = 1, use rainfall data and convert to runoff
% if Bdata.Qf_switch = 2, use riverflow data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default
Bdata.Qf = LochData.Qf * ones(Param.Ndays,1);

if Param.QfData > 0

	% Read rainfall/discharge data
	Qf_data = dlmread(Bdata.Qf_file,',',1,0);
	ncol = size(Qf_data,2);

    % Get latitude and longitude of all data sites
    Date_cols = 4;
    lat = Qf_data(1,Date_cols+1:end);
    lon = Qf_data(2,Date_cols+1:end);

	% Find closest source to current loch
    dist = [];
	for ic = 1:length(lat)
	    dist(ic)=find_dist(LochData.Lat,lat(ic),LochData.Lon,lon(ic),'km');
	end
	icol_Qf = find(dist == min(dist));

	% Select data for required year
	yyyy = Qf_data(:,1);
	iyr = yyyy == Param.Year;
	Bdata.Qf = Qf_data(iyr,Date_cols+icol_Qf);

	% Add spin-up period to data
	Bdata.Qf = [repmat(Bdata.Qf(1),ndays_spin_up,1); Bdata.Qf];

	if Param.QfData == 1
	    % Convert from (mm/day) to river flow (m3/s) assuming conversion 
	    % of 70% of rainfall into river discharge.
	    rainfall = Bdata.Qf;
	    Bdata.Qf = 0.7 * rainfall * LochData.Wshed / (1000 * 86400);
	    % Apply a 3-day smoothing filter
	    a = 1;
	    b = ones(3,1)/3;
	    qff = filter(b,a,Bdata.Qf);
	    Bdata.Qf = qff;
	    % Add direct rainfall onto the water surface
	    Bdata.Qf = Bdata.Qf + rainfall * Hypso.A(1) / (1000 * 86400);
	end

	% Constrain excessively low flows to prevent model instability
	Bdata.Qf(Bdata.Qf < 1) = 1;

    % If data record is shorter than simulation length, extend record using
    % final data value.
    if length(Bdata.Qf) < Param.Ndays
        Bdata.Qf(end:Param.Ndays) = Bdata.Qf(end);
    end

	% Read and print header line
	fid = fopen(Bdata.Qf_file);
	format = repmat('%s',1,ncol);
	header = textscan(fid,format,1,'delimiter',',');
	Source = header{1,Date_cols+icol_Qf};
	fclose(fid);
	disp(['Freshwater discharge data from ',Source{1,1}]);
else
    disp(['Default freshwater discharge = ',num2str(LochData.Qf), ...
        ' cumecs']);
end

% Print river tracer loads
if Param.Ntracer > 0
    C1tmp = LochData.QC1 * 365 * 86400 / 1000;
    disp(['Riverine loading (C1) = ',num2str(C1tmp,'%10.2f'),' tonnes/yr']);
end
if Param.Ntracer > 1
    C2tmp = LochData.QC2 * 365 * 86400 / 1000;
    disp(['Riverine loading (C2) = ',num2str(C2tmp,'%10.2f'),' tonnes/yr']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read wind speed input data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set default
Bdata.Uw = zeros(Param.Ndays,1);

if Param.WindData > 0

    % Read wind speed data
	WS_data = dlmread(Bdata.Wind_file,',',1,0);
	ncol = size(WS_data,2);

	% Get latitude and longitude of all rainfall data sites
	Date_cols = 4;
	lat = WS_data(1,Date_cols+1:end);
	lon = WS_data(2,Date_cols+1:end);

	% Find closest source to current loch
    dist = [];
	for ic = 1:length(lat)
	    dist(ic)=find_dist(LochData.Lat,lat(ic),LochData.Lon,lon(ic),'km');
	end
	icol_WS = find(dist == min(dist));
	% Select data for required year
	yyyy = WS_data(:,1);
	iyr = yyyy == Param.Year;
	Bdata.Uw = WS_data(iyr,Date_cols+icol_WS);
	Ndays = length(Bdata.Uw);

	% Read header line
	fid = fopen(Bdata.Wind_file);
	format = repmat('%s',1,ncol);
	header = textscan(fid,format,1,'delimiter',',');
	Source = header{1,Date_cols+icol_WS};
	fclose(fid);
	disp(['Wind speed source data from ',Source{1,1}]);

	% Add spin-up period to data
	Bdata.Uw = [repmat(Bdata.Uw(1),ndays_spin_up,1); Bdata.Uw];
    
    % If data record is shorter than simulation length, extend record using
    % final data value.
    if length(Bdata.Uw) < Param.Ndays
        Bdata.Uw(end:Param.Ndays) = Bdata.Uw(end);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify coastal salinity and density profiles
% Data from data file, interpolated to 1m
% depth intervals.

% Define the expansion coefficients
Const.alpha = 0.15;
% saline contraction coefficient
Const.beta = 0.78;
% Reference density
Const.rho0 = 1027;

if Param.OBdata > 0
	disp(['Loading ocean boundary data from ',Bdata.OB_file]);
	format = '%n %n %n %n %n %n %n %n%*[^\n]';
	[lon,lat,year,month,day,depth,T,S] = ...
	    textread(Bdata.OB_file,format, ...
             	'delimiter',',','headerlines',1);
	% If the climatology is being used, set the year to Param.year
	id = find(year == Param.Year);
    if year == 0
        year = year + Param.Year;
    elseif length(id) == 0
        disp('ERROR: Mis-match between open boundary data and specified simulation year')
        return
    end

	% Find the closest data position to the specified loch
    dist = [];
	lat2 = lat(depth == 0 & month == 1);
	lon2 = lon(depth == 0 & month == 1);
	for ic = 1:length(lat2)
	    dist(ic)=find_dist(LochData.Lat,lat2(ic),LochData.Lon,lon2(ic),'km');
	end
	ctdlat = lat2(dist == min(dist));
	ctdlon = lon2(dist == min(dist));
	% if there are multiple profiles, choose the first
	ctdlat = ctdlat(1);
	ctdlon = ctdlon(1);
    disp(['Boundary data taken at: ',num2str(ctdlat),' N ',num2str(ctdlon), ...
          ' E']);

	% Get all data for chosen climatology station
	id = lon == ctdlon & lat == ctdlat;
	ctd_z = depth(id);
	ctd_year = year(id);
	ctd_mon = month(id);
	ctd_day = day(id);
	ctd_t = T(id);
	ctd_s = S(id);
	
	% get the dates of each data profile relative to 1st January
	dayno = [datenum(ctd_year,ctd_mon,ctd_day) - datenum(Param.Year,1,1) + 1]';
	iday = dayno(1);
	iprof = 0;

	% For each survey, specify profile of salinity and density
	% and increase vertical spatial resolution to 1-m interval
	% Get profile data
	while iday <= dayno(end)
	    iprof = iprof + 1;
	    id = find(dayno == iday);
	    sal = ctd_s(id);
	    tem = ctd_t(id);
	    z = ctd_z(id);
	    prof_day(iprof) = iday;
	    % Interpolate profiles to 1m depth intervals
	    z_1m = [1:z(end)]';  n_z = length(z_1m);
	    sal_1m(1:n_z,iprof) = interp1(z,sal,z_1m,'cubic');
	    tem_1m(1:n_z,iprof) = interp1(z,tem,z_1m,'cubic');
	    if max(z_1m) < 200
        	sal_1m(n_z+1:200,iprof) = max(sal_1m(n_z,iprof));
        	tem_1m(n_z+1:200,iprof) = max(tem_1m(n_z,iprof));
        	z_1m = [1:200]';
    	end
    	% get next profile
    	if iday < dayno(end); 
	        iday = dayno(id(end)+1); 
	    else
        	iday = dayno(end) + 1;
    	end
	end

	% Increase resolution from monthly to daily intervals. First, if the data 
	% only extend from January - December, extend the
	% boundary data by wrapping around the first and last profiles i.e. insert
	% the December profile at the beginning and add the January profile to the
	% end.
	if ctd_mon(end) == 12
	    prof_day = [prof_day 381];
	    sal_1m = [sal_1m sal_1m(:,1)];
	    tem_1m = [tem_1m tem_1m(:,1)];
	end
	if ctd_mon(1) == 1
	    prof_day = [-15 prof_day];
	    sal_1m = [sal_1m(:,12) sal_1m];
	    tem_1m = [tem_1m(:,12) tem_1m];
	end

	% Interpolate monthly profiles into daily datasets
	dayint = [1:Ndays];
	Bdata.S_ext = interp2(prof_day,z_1m,sal_1m,dayint,z_1m,'linear');
	Bdata.T_ext = interp2(prof_day,z_1m,tem_1m,dayint,z_1m,'linear');

	% Add spin-up period to start of climatology
	Bdata.S_ext = [repmat(Bdata.S_ext(:,1),1,ndays_spin_up) Bdata.S_ext];
	Bdata.T_ext = [repmat(Bdata.T_ext(:,1),1,ndays_spin_up) Bdata.T_ext];

else
    % load boundary data from a matlab file. The file should contain the
    % following variables: Bdata.T0(:,2), Bdata.S0(:,2), Bdata.rho0(:,2).
    load(Bdata.OB_file);    
    S0 = S0'; T0 = T0';
    maxdepth = LochData.Hmax;
    Bdata.S_ext = zeros(maxdepth,size(S0,2));
    Bdata.T_ext = zeros(maxdepth,size(S0,2));
    % Build a layered profile from the layer parameter values
    for ic = 1:size(S0,2);
        ndep1 = min([floor(H0(ic,1)) maxdepth]);
        Bdata.S_ext(1:ndep1,ic) = repmat(S0(1,ic),ndep1,1);
        Bdata.T_ext(1:ndep1,ic) = repmat(T0(1,ic),ndep1,1);
        ndep2 = min([floor(H0(ic,2)) maxdepth]);
        Bdata.S_ext(ndep1+1:ndep2,ic) = repmat(S0(2,ic),ndep2-ndep1,1);
        Bdata.T_ext(ndep1+1:ndep2,ic) = repmat(T0(2,ic),ndep2-ndep1,1);
        ndep3 = max([maxdepth - ndep2 0]);
        Bdata.S_ext(ndep2+1:maxdepth,ic) = repmat(S0(3,ic),ndep3,1);
        Bdata.T_ext(ndep2+1:maxdepth,ic) = repmat(T0(3,ic),ndep3,1);
    end
    % add spin-up period data
    Bdata.S_ext = [repmat(Bdata.S_ext(:,1),1,ndays_spin_up) Bdata.S_ext];
    Bdata.T_ext = [repmat(Bdata.T_ext(:,1),1,ndays_spin_up) Bdata.T_ext];
end

% Derive the external density profile
Bdata.rho_ext = Const.rho0 - Const.alpha * (Bdata.T_ext - 10) ...
    + Const.beta * (Bdata.S_ext - 35);
   
% Specify external tracer concentrations
if Param.Ntracer > 0
    Bdata.C10 = zeros(size(Bdata.S0));
end
if Param.Ntracer > 1
    Bdata.C20 = zeros(size(Bdata.S0));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Surface Heat Flux data if specified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Param.SHFlux == 1
    % Read data
    disp(['Loading surface heat flux data from ',Bdata.SHF_file]);
    SHF_data = dlmread(Bdata.SHF_file,',',1,0);
    ncol = size(SHF_data,2);
    % Select data for required year
    yyyy = SHF_data(:,1);
    iyr = yyyy == Param.Year;
    Bdata.Mslp = SHF_data(iyr,6);
    Bdata.Tair = SHF_data(iyr,7);
    Bdata.Tdew = SHF_data(iyr,8);
    Bdata.Twet = SHF_data(iyr,9);
    Bdata.Relh = SHF_data(iyr,10);
    Bdata.Clcover = SHF_data(iyr,11);
    Bdata.QHin = SHF_data(iyr,12);

    % Add spin-up period to data
    Bdata.Mslp = [repmat(Bdata.Mslp(1),ndays_spin_up,1); Bdata.Mslp];
    Bdata.Tair = [repmat(Bdata.Tair(1),ndays_spin_up,1); Bdata.Tair];
    Bdata.Tdew = [repmat(Bdata.Tdew(1),ndays_spin_up,1); Bdata.Tdew];
    Bdata.Twet = [repmat(Bdata.Twet(1),ndays_spin_up,1); Bdata.Twet];
    Bdata.Relh = [repmat(Bdata.Relh(1),ndays_spin_up,1); Bdata.Relh];
    Bdata.Clcover = [repmat(Bdata.Clcover(1),ndays_spin_up,1); Bdata.Clcover];
    Bdata.QHin = [repmat(Bdata.QHin(1),ndays_spin_up,1); Bdata.QHin];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify variability of coastal density. These variables
% are depth-dependent.
Bdata.sigmarho(1:200,1:Ndays) = zeros;
Bdata.deltaM(1:200,1:Ndays) = zeros;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify the contribution of the M2 + S2 tide to
% the total tidal forcing.
phi = 1.0;
% Specify daily tidal amplitude. Assume a spring-neap 
% cycle where the neap range = 40% of the spring range.
t = [1:Param.Ndays]';
Rmean = 0.7 * LochData.Range;
Ramp = 0.3 * LochData.Range;
Bdata.a0 = phi * 0.5 * (Rmean + Ramp * cos(2 * pi * t / 15));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('ABC_tmp','Bdata', 'Param','Const'); % !da
end
