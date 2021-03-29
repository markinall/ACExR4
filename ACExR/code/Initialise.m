function Initialise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Initialise
%
% Initialises model parameters from topographic and
% boundary forcing data.
%
% Usage:    LochData contains input parameters for loch
%           SillData contains parameters for each sill
%
% Phil Gillibrand
% 14/3/2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('Initialising variables');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Free ("tunable") parameters
% Constant in entrainment velocity calculations
Const.ce = 4.0e-9;
% Intermediary circulation constant (FjordEnv)
% Const.gamma = 17e-4;
% Parameters a and b in the vertical diffusion scheme
Const.gamma = 0.2;
Const.lambda = -2;
% Time step (s)
Const.deltaT = 0.01 * 86400 * Hypso.dz;
% Relaxation time scale (days) for S(2). Currently not used.
Const.relaxT = 0 / (10 * 24 * 3600);
% Gravitational estuarine circulation constant
Const.est = 0.45e5;
% Drag coefficient
Const.Cd = 0.0025;
% Constant in non-stratified flow diffusivity equation
Const.mu = 2.5;
% Constraints on H1
%Const.H1min = max(0.1 * LochData.Hmax, Hypso.z(2));
Const.H1min = Hypso.z(2);
Const.V1min = getV(Const.H1min);
if LochData.Nsill > 0
    Const.H1max = min([0.5*LochData.Hmax SillData.Hmax(1) 20]);
    Const.H1max = min([0.5*LochData.Hmax 20]);
    Const.V1max = getV(Const.H1max);
else
    Const.H1max = 0.5*LochData.Hmax;
    Const.V1max = getV(Const.H1max);
end
% Constraints on H3
if LochData.Nsill > 0
    Const.H3min = LochData.Hmax / 3;
    Const.H3max = LochData.Hmax - SillData.Hmax(1);
    Const.V3max = Hypso.vol(end) - getV(SillData.Hmax(1));
end

% Set parameter values for Etive, 2000 and Creran, 1978 simulations
% obtained from ensemble simulations (July, 2011).
if length(LochData.Name) == 5 & lower(LochData.Name) == 'etive'
    Const.ce = 3.2e-9;
    Const.est = 0.5e5;
    Const.gamma = 0.28;
    Const.lambda = -2.0;
    Const.mu = 1.94;
elseif length(LochData.Name) == 6 & lower(LochData.Name) == 'creran'
    Const.ce = 4.7e-9;
    Const.est = 0.4e5;
    Const.gamma = 0.1;
    Const.lambda = -2.0;
    Const.mu = 2.90;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other constants
% Gravitational acceleration
Const.g = 9.81;
% Von Karmen constant
Const.kappa = 0.4;
% Tidal (M2) period (s)
Const.Tperiod = 44712;
% Day length (seconds)
Const.Tday = 86400;
% Heat capacity of seawater
Const.cp = 4200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error catching
Param.Error = 0;
% Ramp (days)
Const.ramp = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define N and Af
Param.N = Const.ce * (Bdata.Uw .^ 3);
Param.Af = Hypso.A(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Arrays
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise layer parameters.
% Initial conditions depend on whether sills are present.
% If the initial conditions have already been set, do not reset
if Param.ICswitch == 0
    % define parameter arrays
    ndays = Param.Ndays;
    Param.H = zeros(ndays,3);
    Param.V = zeros(ndays,3);
    Param.T = zeros(ndays,3);
    Param.S = zeros(ndays,3);
    Param.rho = zeros(ndays,3);
    Param.Kz = zeros(ndays,2);
    Param.Kzb = zeros(ndays,1);
    Param.DWR = zeros(ndays,1);
    Param.Ri = zeros(ndays,2);
    
    % set parameter values if estuary system contains sills
    if LochData.Nsill == 0
        
        % Set initial estimate based on estuary volume
        Param.V(1,1) = 0.5 * Hypso.vol(end);
        
        % Derive surface and intermediate layer thicknesses
        % Surface layer
        Param.H(1,1) = getH(Param.V(1,1));
        % Check H1 lies in allowable range and adjust volumes if not
        if Param.H(1,1) < Const.H1min
            Param.H(1,1) = Const.H1min;
        elseif Param.H(1,1) > Const.H1max
            Param.H(1,1) = Const.H1max;
        end
        % Re-calculate surface layer volume
        Param.V(1,1) = getV(Param.H(1,1));
        
        % Intermediate layer
        Param.V(1,2) = Hypso.vol(end) - Param.V(1,1);
        vol = Param.V(1,1) + Param.V(1,2);
        Param.H(1,2) = getH(vol) - Param.H(1,1);
        
        % Deep layer thickness is zero with no sills
        Param.H(1,3) = 0;
        Param.V(1,3) = 0;
        
        % Set initial temperature and salinity values
        iz1 = ceil(Param.H(1,1));
        iz2 = min([floor(Param.H(1,1) + Param.H(1,2)), size(Bdata.S_ext,1)]);
        Param.S(1,1) = 0.99 * mean(Bdata.S_ext(1:iz1,1));
        Param.S(1,2) = mean(Bdata.S_ext(iz1:iz2,1));
        Param.S(1,3) = 0;
        Param.T(1,1) = mean(Bdata.T_ext(1:iz1,1));
        Param.T(1,2) = mean(Bdata.T_ext(iz1:iz2,1));
        Param.T(1,3) = 0;
        if Param.Ntracer > 0
            Param.C1(1,1) = mean(Bdata.C1_ext(1:iz1,1));
            Param.C1(1,2) = mean(Bdata.C1_ext(iz1:iz2,1));
            Param.C1(1,3) = 0;
        end
        if Param.Ntracer > 1
            Param.C2(1,1) = mean(Bdata.C2_ext(1:iz1,1));
            Param.C2(1,2) = mean(Bdata.C2_ext(iz1:iz2,1));
            Param.C2(1,3) = 0;
        end
    else
        % With at least one sill
        % Set layer thicknesses
        Param.H(1,1) = 0.5 * SillData.Hmax(1);
        Param.H(1,2) = 0.5 * SillData.Hmax(1);
        Param.H(1,3) = LochData.Hmax - SillData.Hmax(1);
        Param.V(1,1) = getV(Param.H(1,1));
        h2 = min([Param.H(1,1)+Param.H(1,2), Hypso.z(end)]);
        Param.V(1,2) = getV(h2) - Param.V(1,1);
        Param.V(1,3) = Hypso.vol(end) - Param.V(1,1) - Param.V(1,2);
        iz1 = ceil(Param.H(1,1));
        iz2 = floor(Param.H(1,1) + Param.H(1,2));
        iz3 = floor(Param.H(1,1) + Param.H(1,2) + Param.H(1,3));
        % Salinity
        Param.S(1,1) = 0.99 * mean(Bdata.S_ext(1:iz1,1));
        Param.S(1,2) = mean(Bdata.S_ext(iz1:iz2,1));
        Param.S(1,3) = mean(Bdata.S_ext(iz2:iz3,1));
        % Temperature
        Param.T(1,1) = mean(Bdata.T_ext(1:iz1,1));
        Param.T(1,2) = mean(Bdata.T_ext(iz1:iz2,1));
        Param.T(1,3) = mean(Bdata.T_ext(iz2:iz3,1));
        % Tracer 1
        if Param.Ntracer > 0
            Param.C1(1,1) = mean(Bdata.C1_ext(1:iz1,1));
            Param.C1(1,2) = mean(Bdata.C1_ext(iz1:iz2,1));
            Param.C1(1,3) = mean(Bdata.C1_ext(iz2:iz3,1));
        end
        % Tracer 2
        if Param.Ntracer > 1
            Param.C2(1,1) = mean(Bdata.C2_ext(1:iz1,1));
            Param.C2(1,2) = mean(Bdata.C2_ext(iz1:iz2,1));
            Param.C2(1,3) = mean(Bdata.C2_ext(iz2:iz3,1));
        end
    end
else
    % Initialise parameter arrays
    % T, S, H already set
    ndays = Param.Ndays;
    Param.H(2:ndays,:) = zeros(ndays-1,3);
    Param.T(2:ndays,:) = zeros(ndays-1,3);
    Param.S(2:ndays,:) = zeros(ndays-1,3);
    Param.V = zeros(ndays,3);
    Param.Kz = zeros(ndays,2);
    Param.Kzb = zeros(ndays,1);    
    Param.DWR = zeros(ndays,1);
    Param.Ri = zeros(ndays,2);
end
Param.rho = Const.rho0 - Const.alpha * (Param.T - 10) + ...
    Const.beta * (Param.S - 35);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = getH(V)

global Hypso

iz1 = find(Hypso.vol < V); iz1 = iz1(end);
iz2 = find(Hypso.vol >= V); iz2 = iz2(1);
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
iz2 = ceil(H / Hypso.dz) + 1;
f = (H - Hypso.z(iz1)) / Hypso.dz;
A = (1-f)*Hypso.A(iz1) + f*Hypso.A(iz2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
