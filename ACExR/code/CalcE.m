function CalcE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function CalcE
%
% The main calculation module which determines the layer depths
%   and volumes and their exchange rates.
%
% Usage:    LochData contains topographic and other data
%           Bdata contains boundary forcing data
%           Param contains various parameters (cf Initialise.m)
%           E contains the main output exchange rate data for
%               each layer.
%
%           The internal structured array D contains the most recent values
%           of all temporal variables (i.e. layer parameters and fluxes),
%           calculated within the sub-daily time step. The final sub-daily
%           values are passed to Param, and the sub-daily fluxes are
%           integrated to give daily-average values in E.
%
%           The final 
% Phil Gillibrand
% 27/3/2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('Calculating exchange rates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get length of simulation
Ndays = Param.Ndays;

% Number of time steps per day
Nsteps_per_day = 86400 / Const.deltaT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise arrays
E.Qt = zeros(Ndays,1);
E.Qe = zeros(Ndays,1);
E.Qh = zeros(Ndays,2);
E.Qi = zeros(Ndays,1);
E.Qw12 = zeros(Ndays,1);
E.Qw23 = zeros(Ndays,1);
E.Kz12 = zeros(Ndays,1);
E.Kz23 = zeros(Ndays,1);
E.Kzb = zeros(Ndays,2);
E.Ut = zeros(Ndays,1);
E.Ue = zeros(Ndays,2);
Bdata.Qshf = zeros(Ndays,1);
% Initialise variables
D.Qe = 0;
D.Qh = [0 0];
D.Qi = 0;
D.Qw12 = 0;
D.Qw23 = 0;
D.Kz12 = 0;
D.Kz23 = 0;
D.Kzb = 0;
D.Qshf = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display constraints on H and V
disp(['H1 min and max: ',num2str([Const.H1min Const.H1max]),'  m'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign the freshwater flux to the exchange rate array
E.Qf = Bdata.Qf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise layer parameters for Day 1.
day = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate surface and intermediate layer thicknesses
D.H = Param.H(1,:);
if LochData.Nsill > 0
    if D.H(3) ~= LochData.Hmax - D.H(1) - D.H(2)
        D.H(3) = LochData.Hmax - D.H(1) - D.H(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign surface layer temperature and salinity and density based on 
% external boundary conditions and mixing/buoyancy balance in surface 
% layer.
D.T = Param.T(1,:);
D.S = Param.S(1,:);
if Param.Ntracer > 0
    D.C1 = Param.C1(1,:);
end
if Param.Ntracer > 1
    D.C2 = Param.C2(1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify density
alpha = Const.alpha;
beta = Const.beta;
D.rho = Const.rho0 - alpha * (D.T - 10) + beta * (D.S - 35);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate tidal exchange efficiency
eps = Efficiency;
Param.eps = eps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise time
day0 = 0;
Param.t = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If hot-starting, load variables
if Param.hotstart == 1
    disp(['Hot-start: loading data from ',Param.hotstartfile])
    load(Param.hotstartfile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ')
disp(['Initial Layer thicknesses : ',num2str(D.H)])
disp(['Initial Layer salinities : ',num2str(D.S)])
disp(['Initial Layer temperatures : ',num2str(D.T)])
disp(['Initial Layer densities : ',num2str(D.rho)])
disp(' ')
disp(['Run length : ',num2str(Ndays),' days'])
disp(' ')
%disp('N.B. Entrainment is switched OFF');
%disp(' ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Loop.
for day = 1:Ndays

    if rem(day,10) == 0; disp(['Day ',num2str(day)]); end
    Param.day = day + day0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate daily exchange rates and estimate new mixed layer
    % depth based on daily wind mixing and buoyancy input.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set open boundary temperature and salinities
    iz1 = ceil(D.H(1));
    if LochData.Nsill > 0
        iz2 = floor(SillData.Hmax(1));
    else
        iz2 = min([floor(D.H(1,1) + D.H(1,2)) size(Bdata.S_ext,1)]);
    end
    Bdata.T0(day,1) = mean(Bdata.T_ext(1:iz1,day));
    Bdata.S0(day,1) = mean(Bdata.S_ext(1:iz1,day));
    Bdata.T0(day,2) = Bdata.T_ext(iz2,day);
    Bdata.S0(day,2) = Bdata.S_ext(iz2,day);
    if Param.Ntracer > 0
        Bdata.C10(day,1) = mean(Bdata.C1_ext(1:iz1,day));
        Bdata.C10(day,2) = Bdata.C1_ext(iz2,day);
    end
    if Param.Ntracer > 1
        Bdata.C20(day,1) = mean(Bdata.C2_ext(1:iz1,day));
        Bdata.C20(day,2) = Bdata.C2_ext(iz2,day);
    end
    T0 = Bdata.T0(day,:);
    S0 = Bdata.S0(day,:);
    rho0 = Const.rho0 - alpha * (T0 - 10) + beta * (S0 - 35);
    Bdata.rho0(day,:) = rho0;
    if Param.Ntracer > 0
        C10 = Bdata.C10(day,:);
    end
    if Param.Ntracer > 1
        C20 = Bdata.C20(day,:);
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate tidal exchange rate.
    E.Qt(day,1) = 2 * eps * Bdata.a0(day) * Param.Af / Const.Tperiod;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the exchange rate due to the intermediary circulation
    if LochData.Nsill == 0;
        iz1 = max(round(LochData.Hmean), 1);
        deltaM = Bdata.deltaM(iz1);
        E.Qi(day,1) = Const.gamma * sqrt(Const.g * Hypso.vol(end) ...
                    * Param.Af * deltaM / (LochData.Len * D.rho(2)));
    else
        iz1 = max(round(SillData.Hmax(1)), 1);
        deltaM = Bdata.deltaM(iz1);
        E.Qi(day,1) = Const.gamma * sqrt(Const.g * SillData.Xarea(1) ...
                    * Param.Af * deltaM / D.rho(2));
        if Param.Af / SillData.Xarea(1) > 1e4
            E.Qi(day,1) = (1/6) * SillData.Xarea(1) * sqrt(2 * Const.g ...
                        * deltaM / D.rho(2));
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Estimate background mixing due to internal wave activity
    if LochData.Nsill > 0
        W = BasinMixing;
        % W_sum = W(1); % commented out 8 Nov 2011
        W_sum = sum(W); % inserted on 8 Nov 2011
        if length(W) > 1; W_sum = W_sum + sum(W(2:end)); end
        delrho = max(D.rho(3) - D.rho(1), 0.1);
        Kz_back = W_sum / (Const.g * delrho);
        %Kz_back = W_sum / Const.g;
    else
        Kz_back = 0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step forward in time, calculating layer parameters
    for istep = 1:Nsteps_per_day

        % time 
        Param.t = Param.t + Const.deltaT;
        ramp = min(Param.t/(Const.ramp * 86400), 1);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Ensure layer volumes and cross-secional and interfacial areas
        % between layers are consistent with layer thicknesses.
        % Surface layer
        D.V(1) = getV(D.H(1));
        D.xarea(1) = getX(D.H(1));
        D.A_ifc(1) = getA(D.H(1));
        % Intermediate layer
        h2 = min([D.H(1)+D.H(2), Hypso.z(end)]);
        D.V(2) = getV(h2) - D.V(1);
        D.xarea(2) = getX(h2) - D.xarea(1);
        D.A_ifc(2) = getA(h2);
        % Deep layer
        if LochData.Nsill > 0
            D.V(3) = Hypso.vol(end) - D.V(1) - D.V(2);
        else
            D.V(3) = 0;
        end

        % Reset maximum permitted surface layer thickness and volume
        Const.H1max = 0.5 * (LochData.Hmax - D.H(3));
        Const.V1max = getV(Const.H1max);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Look to see whether external density is greater than deep water
        % density. If it is, the deep water will be replaced. Switch on
        % deep water renewal switch.
        D.DWR = 0;
        if LochData.Nsill > 0
            if rho0(2) > D.rho(3) & D.H(3) < Const.H3max;
                D.DWR = 1;
            end
        end
        Param.DWR(day,1) = Param.DWR(day) + D.DWR * Const.deltaT / 86400;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get new value of M
        Param.M = Const.g * Const.beta * D.S(1) / Const.rho0;

        % Calculate the heat flux through the water surface (W/m2)
        if Param.SHFlux == 1; 
            D.Qshf = SurfaceHeatFlux; 
        else
            D.Qshf = 0;
        end

        % get entrainment due to wind stirring
        D.Qh = ramp * WindStirring;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate gravitational estuarine circulation, defined as a function
        % of the horizontal salinity gradient. See Li et al. (1999).
        D.Qe = ramp * Const.est * (rho0(1) - D.rho(1)) / Const.rho0;
        if D.Qe < 0; D.Qe = 0; end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Diagnose turbulent diffusivity and entrainment between layers 1 and 2
        KZ = Kz12_general;
        %delrho = max(D.rho(2) - D.rho(1), 0.1);
        Param.Kzb(day,1) = Param.Kzb(day,1) + Kz_back * Const.deltaT ...
                / 86400;
        D.Kzb(1) = ramp * Kz_back * D.A_ifc(1) / sum(D.H);
        D.Kz12 = ramp * KZ;

        % Calculate tidal entrainment.
        % N.B. No entrainment is modelled during deep water 
        % renewal events.
        Ent = Entrainment;
        %Ent = [0 0];
        D.Qw12 = ramp * Ent(1);
        D.Qw23 = ramp * Ent(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Diagnose turbulent diffusivity and entrainment between layers 
        % 2 and 3 if layer 3 is present.
        if LochData.Nsill > 0
            KZ = Kz23_general;
            % Add in background mixing due to internal wave activity
            D.Kzb(2) = ramp * Kz_back * D.A_ifc(2) / sum(D.H);
            D.Kz23 = ramp * KZ + D.Kzb(2);
        end        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Record present layer thicknesses and volumes
        D.Hp = D.H;
        D.Vp = D.V;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the new layer thicknesses and volumes
        V_new(1) = D.V(1) + Const.deltaT * (D.Qh(1) + D.Qw12);
        V_new(2) = D.V(2) - Const.deltaT * (D.Qh(1) - D.Qh(2) + D.Qw12 ...
                                          - D.Qw23);
        V_new(3) = D.V(3) - Const.deltaT * (D.Qh(2) + D.Qw23);
        if LochData.Nsill > 0 & V_new(3) > Const.V3max
            dV = V_new(3) - Const.V3max;
            V_new(3) = Const.V3max;
            V_new(2) = V_new(2)+ dV;
            D.Qe =(Const.V3max - D.V(3)) / Const.deltaT;
        end
                                              
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Derive surface and intermediate layer thicknesses
        % Surface layer
        H_new(1) = getH(V_new(1));
        % Intermediate layer
        vol = min([V_new(1) + V_new(2), Hypso.vol(end)]);
        H_new(2) = getH(vol) - H_new(1);
        % Deep layer
        if LochData.Nsill > 0
            H_new(3) = LochData.Hmax - H_new(1) - H_new(2);
        else
            H_new(3) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update layer thickness and volume parameters
        D.H = H_new;
        D.V = V_new;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the x-sectional and interfacial areas for layers one and two
        % Get the x-sectional area of layer one
        D.xarea(1) = getX(D.H(1));
        % Get the interfacial area between layers one and two
        D.A_ifc(1) = getA(D.H(1));
        % Get the x-sectional area of layer two
        h2 = min([D.H(1)+D.H(2), Hypso.z(end)]);
        D.xarea(2) = getX(h2) - D.xarea(1);
        % Get the interfacial area between layers two and three
        D.A_ifc(2) = getA(h2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate layer temperatures and salinities
        Tnew = TracerSI(D.T,D.T(1),T0);
        if Param.SHFlux == 1; 
            Hbar = 0.5 * (D.Hp(1) + D.H(1));
            Tnew(1) = Tnew(1) + Const.deltaT * D.Qshf / ...
                (Const.cp * D.rho(1) * Hbar);
        end
        [Snew] = TracerSI(D.S,0,S0);
        
        % Attempt to catch errors
        if Tnew(1) < 0 | Tnew(1) > 30
            Param.Error = 1;
            return;
        end
        if Tnew(2) < 0 | Tnew(2) > 30
            Param.Error = 2;
            return;
        end
        if Snew(1) < 0 || Snew(1) > 40
            Param.Error = 3;
            return;
        end
        if Snew(2) < 0 | Snew(2) > 40
            Param.Error = 4;
            return;
        end
     
        % Nutrients.
        % Convert nutrient mass flux (kg/s) to concentration (kg/m3)
        if Param.Ntracer > 0 & day > Param.ndays_spin_up
            C1river = LochData.QC1 / Bdata.Qf(day);
            [C1new] = TracerSI(D.C1,C1river,C10);
            D.C1 = C1new;
        end
        if Param.Ntracer > 1 & day > Param.ndays_spin_up
            C2river = LochData.QC2 / Bdata.Qf(day);
            [C2new] = TracerSI(D.C2,C2river,C20);
            D.C2 = C2new;
        end
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update all layer variables
        D.T = Tnew;
        D.S = Snew;
        % Derive density
        D.rho = Const.rho0 - alpha * (D.T - 10) + beta * (D.S - 35);
                            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Integrate daily fluxes
        E.Qe(day) = E.Qe(day) + D.Qe;
        E.Qh(day,:) = E.Qh(day,:) + D.Qh;
        E.Qw12(day) = E.Qw12(day) + D.Qw12;
        E.Qw23(day) = E.Qw23(day) + D.Qw23;
        E.Kz12(day) = E.Kz12(day) + D.Kz12;
        E.Kz23(day) = E.Kz23(day) + D.Kz23;
        E.Kzb(day,:) = E.Kzb(day,:) + D.Kzb;
        E.Ue(day,1) = E.Ue(day,1) + (D.Qe + E.Qi(day)) / D.xarea(1);
        E.Ue(day,2) = E.Ue(day,2) + (D.Qe + E.Qi(day)) / D.xarea(2);
        E.Ut(day,1) = E.Ut(day,1) + E.Qt(day) / D.xarea(2);
        Bdata.Qshf(day) = Bdata.Qshf(day) + D.Qshf;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign end-of-day values to output arrays.
        if istep == 1
            Param.T(day,:) = zeros;
            Param.S(day,:) = zeros;
            Param.H(day,:) = zeros;
            Param.V(day,:) = zeros;
            if Param.Ntracer > 0
                Param.C1(day,:) = zeros;
            end
            if Param.Ntracer > 1
                Param.C2(day,:) = zeros;
            end
        end
        
        Param.H(day,:) = Param.H(day,:) + D.H / Nsteps_per_day;
        Param.V(day,:) = Param.V(day,:) + D.V / Nsteps_per_day;
        Param.T(day,:) = Param.T(day,:) + D.T / Nsteps_per_day;
        Param.S(day,:) = Param.S(day,:) + D.S / Nsteps_per_day;
        if Param.Ntracer > 0; 
             Param.C1(day,:) = Param.C1(day,:) + D.C1 / Nsteps_per_day; 
        end
        if Param.Ntracer > 1; 
             Param.C2(day,:) = Param.C2(day,:) + D.C2 / Nsteps_per_day; 
        end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derive layer densities
    Param.rho(day,:) = Const.rho0 - alpha * (Param.T(day,:) - 10) ...
                        + beta * (Param.S(day,:) - 35);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate daily mean fluxes and velocities
    E.Qe(day) = E.Qe(day) / Nsteps_per_day;
    E.Qh(day,:) = E.Qh(day,:) / Nsteps_per_day;
    E.Qw12(day) = E.Qw12(day) / Nsteps_per_day;
    E.Qw23(day) = E.Qw23(day) / Nsteps_per_day;
    E.Kz12(day) = E.Kz12(day) / Nsteps_per_day;
    E.Kz23(day) = E.Kz23(day) / Nsteps_per_day;
    E.Kzb(day,:) = E.Kzb(day,:) / Nsteps_per_day;
    E.Ue(day,:) = E.Ue(day,:) ./ Nsteps_per_day;
    E.Ut(day) = E.Ut(day) / Nsteps_per_day;
    Bdata.Qshf(day) = Bdata.Qshf(day) / Nsteps_per_day;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save results to file each day
    save results E Param Bdata SillData LochData Const Hypso
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save final set of results to file
Save_results;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save variables for hot-start
save hotstart.mat D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = getH(V)

global Hypso

iz1 = find(Hypso.vol < V); 
if isempty(iz1)
    iz1 = 1;
else
    iz1 = iz1(end);
end
iz2 = find(Hypso.vol >= V);
if isempty(iz2)
    iz2 = length(Hypso.z);
else
    iz2 = iz2(1);
end
f = (V - Hypso.vol(iz1))/(Hypso.vol(iz2) - Hypso.vol(iz1));
H = f * Hypso.z(iz2) + (1-f)*Hypso.z(iz1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function V = getV(H)

global Hypso

iz1 = floor(H / Hypso.dz) + 1;
iz2 = ceil(H / Hypso.dz) + 1;
f = (H - Hypso.z(iz1)) / Hypso.dz;
V = (1-f)*Hypso.vol(iz1) + f*Hypso.vol(iz2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = getX(H)

global Hypso

iz1 = floor(H / Hypso.dz) + 1;
iz2 = ceil(H / Hypso.dz) + 1;
f = (H - Hypso.z(iz1)) / Hypso.dz;
X = (1-f)*Hypso.xarea(iz1) + f*Hypso.xarea(iz2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = getA(H)

global Hypso

iz1 = floor(H / Hypso.dz) + 1;
if iz1 > length(Hypso.z); iz1 = length(Hypso.z); end;
iz2 = ceil(H / Hypso.dz) + 1;
if iz2 > length(Hypso.z); iz2 = length(Hypso.z); end;
f = (H - Hypso.z(iz1)) / Hypso.dz;
A = (1-f)*Hypso.A(iz1) + f*Hypso.A(iz2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
