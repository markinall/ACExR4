function KZ = Kz12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function KZ = Kz12
%
% evaluates definite symbolic integral of the Kz(x) along the axis
% of a basin. Kz is based on a general Richardson-number dependent
% formulation used by Munk & Anderson (1948), Pacanowski & Philander 
% (1981) among others.
%
% called from CalcE function.
%
% for reference:
%
% Kv = Kv0(1 + aRi)^b
%
% Kv0 = Cd * U * H2 * Sc, where Cd is a drag coefficient, U is the 
% velocity in the intermediate layer, H2 is the intermediate layer 
% thickness, and Sc is a Schmidt Number.
%
% a, b are constants e.g.:
% a = 3.33, b = -1.5 (Munk & Anderson, 1948)
% a = 5, b = -3 (Pacanowski & Philander, 1981)
% a = free, b = -1 (analagous to Babson et al., 2006)
%
% Phil Gillibrand
% June 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the constants
a = Const.gamma;
b = Const.lambda;

% Get required parameters
day = Param.day;
% Estuarine exchange
Qe = abs(D.Qe); 
% Tidal exchange
Qt = E.Qt(day);
% Salinity difference between layers 2 and 3
DelS = D.S(2) - D.S(1); 
% Distance between layer centres
DiffZ = 0.5 * (D.H(1) + D.H(2));
% Linear EoS conversion to density
DelSig = Const.beta * DelS; 
% Rough mean density in kgm-3
MeanSig = 1025; 
% Approximate buoyancy frequency squared (s-2)
N_sq = (Const.g/MeanSig) * (DelSig/DiffZ); 
N_sq = min([max([N_sq 1e-6]) 1e-2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the x-sectional area of layer one
xarea_1 = D.xarea(1);
% Get the x-sectional area of layer two
xarea_2 = D.xarea(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive estuarine induced flow difference between surface
% and intermediate layers +ve
U1e = max([Qe / xarea_1 0.01]);
U2e = max([Qe / xarea_2 0.01]); 
% Derive tidally induced flow in intermediate layer +ve
U2t = max([Qt / xarea_2 0.01]);

% Calculate the sum velocity and the velocity shear
Unet = U2e + U2t + U1e;
Ushr = Unet / DiffZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the overall (bulk) Richardson Number
Rio = max([N_sq / (Ushr * Ushr) -4.5]);
Param.Ri(day,1) = Param.Ri(day,1) + Rio * Const.deltaT / 86400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the vertical eddy diffusion coefficient
Kv0 = Const.mu * Const.Cd * Unet * min([D.H(1) D.H(2)]);
KZ = Kv0 * (1 + a * Rio)^b;
Param.Kz(day,1) = Param.Kz(day,1) + KZ * Const.deltaT / 86400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the interfacial area between layers one and two
A_12 = D.A_ifc(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to a flux of m3/s
KZ = KZ * A_12 / DiffZ;

% end function
