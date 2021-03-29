function Q = SurfaceHeatFlux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Q = SurfaceHeatFlux(Ts);
%
% Calculates surface heat flux to surface layer.
%
% Usage: Q (W/m2) is the net heat flux INTO the surface layer      
%        Ts is the current surface layer temperature
%
% Phil Gillibrand
% 7/08/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Day
day = Param.day;

% Surface Temperature
Ts = D.T(1);

% Assign local variables
ur = Bdata.Uw(day);
z = Bdata.StHgt;
Ta = Bdata.Tair(day);
Tw = Bdata.Twet(day);
rh = Bdata.Relh(day) * 100;
Pa = Bdata.Mslp(day);
nc = Bdata.Clcover(day) / 8;

% Incident radiation
QI = Bdata.QHin(day);

% Latent and Sensible heat fluxes
A = hfbulktc(ur,z,Ta,z,rh,z,Pa,Ts);

% Sensible heat flux (W/m2) into the ocean
QS = A(1);

% Latent heat flux (W/m2) into the ocean
QL = A(2);

% Long-wave radiation
sigma = 5.7e-8;
Tskelv = Ts + 273;
ew = 6.112 * exp(17.67 * Tw / (Tw + 243.5));
ea = ew - Pa * (Ta - Tw) * 0.00066 * (1 + 0.00115 * Tw);
QB = 0.985 * sigma * (Tskelv ^ 4) * (0.39 - 0.05 * sqrt(ea)) * ...
                                            (1 - 0.6 * nc * nc);

% Net heat flux (W/m2)
Q = QI + QS + QL - QB;

end



