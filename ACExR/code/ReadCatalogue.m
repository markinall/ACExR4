function ReadCatalogue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ReadCatalogue;
%
% Reads physical characteristics for 'Lochname' from 
% sea loch catalogue database (in csv format).
% All data converted to SI units.
%
% Usage:    Lochname is a string containing the loch name
%
% Phil Gillibrand
% 10/3/2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lochname = LochData.Name;
disp(' ');
disp(['Searching sea loch database for Loch ',Lochname]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name of sea loch catalogue database file
filename = '../Catalogue/Sealochs.csv';
% specify data format
format='%n%s%s%s%n%n%s%s';
for ic = 1:19; format = [format,'%n']; end

% Read data file
[ID,LochNames,Comment,OS_ref,Lat,Lon,Admiralty,Admiralt_1, ...
data(:,1),data(:,2),data(:,3),data(:,4),data(:,5), ...
data(:,6),data(:,7),data(:,8),data(:,9),data(:,10), ...
data(:,11),data(:,12),data(:,13),data(:,14),data(:,15), ...
data(:,16),data(:,17),data(:,18),data(:,19)] = ...
textread(filename,format,'headerlines',1,'delimiter',',');    

LochNamesDouble = double(upper(char(LochNames)));
[n m] = size(LochNamesDouble);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the data row for specified 'loch' by converting 
% names to double precision values and matching.
lochdouble = double(upper(Lochname));
nchar = length(lochdouble);
lochdouble(nchar+1:m) = ones * 32;

% Search through database for specified 'loch'
for irow = 1:n
    if lochdouble == LochNamesDouble(irow,:)
        row = irow;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign to structured array 'LochData'
LochData.Lat = Lat(row);
LochData.Lon = Lon(row);
LochData.Len = data(row,1) * 1e3;
LochData.Range = data(row,2);
LochData.Hmax = data(row,3);
LochData.HWarea = data(row,4) * 1e6;
LochData.LWarea = data(row,5) * 1e6;
LochData.Area2m = data(row,6) * 1e6;
LochData.Area5m = data(row,7) * 1e6;
LochData.Area10m = data(row,8) * 1e6;
LochData.LWvol = data(row,9) * 1e6;
LochData.Wshed = data(row,10) * 1e6;
LochData.Rain = data(row,11) * 1e-3;
LochData.Qf = data(row,12) * 1e6 / (365 * 24 * 3600);
LochData.Hmean = data(row,14);
LochData.Nsill = data(row,17);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set river tracer concentrations to zero if required
if Param.Ntracer > 0
    LochData.QC1 = 0;
end
if Param.Ntracer > 1
    LochData.QC2 = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the HW volume
LochData.HWvol = LochData.LWvol + 0.5 * (LochData.HWarea ...
    + LochData.LWarea) * LochData.Range;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate flushing time
LochData.Tf = 89424 * LochData.LWvol / ...
    ((LochData.HWarea + LochData.LWarea)*0.7*LochData.Range);

clear data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LochData.Nsill > 0
    % Get sill data for specified loch
    % First find rows containing correct data
    filename = '../Catalogue/Sills.csv';
    format='%n%n%s%n%n%n%n%n%n%n%n%n%n%*[^\n]';
    [lochid,sillno,c, data(:,1), ...
    data(:,2),data(:,3),data(:,4),data(:,5),data(:,6), ...
    data(:,7),data(:,8),data(:,9),data(:,10)] = ...
    textread(filename,format,'headerlines',1,'delimiter',',');    
    irow = find(lochid == row);
    range = min(irow):max(irow);

    % Assign to structured array 'SillData'
    SillData.Len = data(range,1);
    SillData.HWwid = data(range,2);
    SillData.LWwid = data(range,3);
    SillData.Hmax = data(range,4);
    SillData.Hmean = data(range,5);
    SillData.Xarea = data(range,6);
    SillData.HWuparea = data(range,7) * 1e6;
    SillData.LWuparea = data(range,8) * 1e6;
    SillData.Current = data(range,9) * 1e-2;
    SillData.Hbasin = data(range,10);
    Param.Bm = SillData.Xarea(1) / SillData.Hmean(1);
else
    SillData.Len = 0;
    SillData.HWwid = 0;
    SillData.LWwid = 0;
    SillData.Hmax = 0;
    SillData.Hmean = 0;
    SillData.Xarea = 0;
    SillData.HWuparea = 0;
    SillData.LWuparea = 0;
    SillData.Current = 0;
    SillData.Hbasin = 0;
    Param.Bm = (LochData.HWarea + LochData.LWarea) / (2 * LochData.Len);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modify the variables from the sea loch catalogue such that the 
% datum is Mean Sea Level rather than the Low Water mark.
LochData.Hmean = LochData.Hmean + LochData.Range / 2;
LochData.Hmax = LochData.Hmax + LochData.Range / 2;
if LochData.Nsill > 1
    SillData.Hmax = SillData.Hmax + LochData.Range / 2;
    SillData.Hmean = SillData.Hmean + LochData.Range / 2;
    SillData.Xarea = SillData.Xarea + 0.5 * LochData.Range ...
                   * (SillData.HWwid + 3 * SillData.LWwid) / 4;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get estimated surface areas of individual basins
filename_s = '../Catalogue/sill_basin_areas.csv';
format = '%s %n %n';
[loch, a1, a2] = ...
    textread(filename_s,format,'delimiter',',','headerlines',1);
alllochs = double(upper(char(loch)));
[n, m] = size(alllochs);
thisloch = double(upper(Lochname));
nchar = length(thisloch);
thisloch(nchar+1:m) = ones * 32;
row = [];
for irow = 1:n
    if thisloch == alllochs(irow,:)
        row = [row irow];
    end
end
if isempty(row) == 1
    SillData.Barea = 0;
else
    SillData.Barea = a2(row) * 1e6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end