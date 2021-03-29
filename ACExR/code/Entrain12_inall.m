function Ent12 = Entrain12_inall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ent12 = Entrain12_inall
% numerical evaluation of entrainment from layer 1 into layer 2
% by barotropic tidal current subducted under layer 1
% Based on the Ri parameterisation of Ellison and Turner
%
% called from CalcE function.
%
% for reference: A, B, and C are tunable constants for the Turner scheme.
%
% A=0.08;,B=0.1;,C=5;
% Mark Inall
% 29 March 2006

%%%%%%%%%% MEI Comment 14th Dec 2007 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error corrected in orogonal code - factor of U missing in entrainment
% velocity. the code now runs with the A,B,C values as given by Turner
%correct values from Turner 1986 are: A=0.08;,B=0.1;,C=5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEI Comment - Feb 2008 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entraimnet param of Princevac et al, JFM 2005, vol533 pp259-268
% implemented - Turner 1986 param no longer used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MEI - 12 Feb 2008 comment
% shear between L1 and L2 now includes the esyuarine velocity of the upper
% layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get required parameters
day = Param.day;
% Estuarine exchange
Qe = D.Qe; 
% Tidal exchange
Qt = E.Qt(day);
% Intermediary circulation
Qi = E.Qi(day);
% Vertical distance between layer centres
DiffZ = 0.5*(D.H(1) + D.H(2));
% Square of dz for velocity shear
dzsq=DiffZ^2;
% Density difference between layers 1 and 2.
DelSig = D.rho(2) - D.rho(1);
% Rough mean density in kgm-3
MeanSig = 0.5 * (D.rho(1) + D.rho(2)); 
% Approximate buoyancy frequency squared (s-2)
N_sq = (Const.g/MeanSig) * (DelSig/DiffZ); 
N_sq = max([N_sq 1e-6]);
% Length of loch
L = LochData.Len; 
% Get the x-sectional area of layer one
xarea_1 = D.xarea(1);
% Get the x-sectional area of layer two
xarea_2 = D.xarea(2);
% Get the interfacial area between layers one and two
A_ifc = D.A_ifc(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive estuarine induced flow in surface layer
U1e = (Qe + Qi) / xarea_1; 
% Derive estuarine induced flow in intermediate layer
U2e = (Qe + Qi) / xarea_2; 
% Maximum tidal velocity over the sill
if LochData.Nsill > 0
    % U0t = Qt / SillData.Xarea(1);
    U0t = pi * Qt / (2 * Param.eps * SillData.Xarea(1));
    U0e = 2 * Qe / SillData.Xarea(1);
else
    % U0t = Qt / (Hypso.vol(end) / LochData.Len);
    U0t = pi * Qt / (2 * Param.eps * (xarea_1 + xarea_2));
    U0e = Qe / xarea_2;
end
% Amplitude of tidal velocity in layer 2
U2t = pi * Qt / (2 * Param.eps * xarea_2);
% Length scale of velocity decay in layer 2
Xl = L/5;
%if U0t > 0; Xl = L * U2t / U0t; end
% M2 Tidal period
T=Const.Tperiod;
omega = 2*pi/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate entrainment for the period from the present time step 
% to the next time step.
t=Param.t;
U2 = U2e + U2t*sin(omega*t);
x=0:20:L-10;
w=zeros(1,length(x));
if U2 > 0
    U2x = U2e + U2t*sin(omega*t)*exp(-x/Xl);
    %U2x = U2 * exp(-x/Xl);
    dudz2 = (U2x .* U2x / dzsq);
    Ri = N_sq ./ dudz2;
    % Turner (1986) formulation
    %id = Ri < 0.8;
    %w(id) = U2x(id) .* (A - B * Ri(id)) ./ (1 + C * Ri(id));
    % Princevac et al (2005) formulation
     Ri = max(Ri, 0.15);
     id = Ri < 1.5;
     w(id) = U2x(id) .* (0.05 * Ri(id).^(-0.75));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the entrainment to prevent loss of surface layer volume
% when surface layer is already at its minimum value.
Q_scale = (D.H(1) - Const.H1min) / (Const.H1max - Const.H1min);
Q_scale = max([min([Q_scale 1]) 0]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entrainment velocity from layer 1 into layer 2
w12 = mean(w);
% entrainment flux in m3s-1
% N.B. Vertical fluxes positive upwards
Ent12 = -A_ifc * Q_scale * w12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
