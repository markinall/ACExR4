function Ent23 = Entrain23_inall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ent23 = Entrain23_inall
% numerical evaluation of entrainment from layer 3 into layer 2
% by barotropic tidal current.
% Based on the Ri parameterisation of Ellison and Turner
%
% called from CalcE function.
%
% for reference: A, B, and C are tunable constants for the Turner scheme.
%
% A=0.00008; B=0.00; C=10;
% A=0.08; B=0.1; C=5; % From Turner 1986

% Entrain23 by Phil Gillibrand
% November 2006
% Based on Entrain12 code by Mark Inall
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get required parameters
day = Param.day;
% Maximum thickness of layer 3
H3_max = min([Const.H3max LochData.Hmax-2*D.H(1)]);
% Estuarine exchange
Qe = D.Qe; 
% Tidal exchange
Qt = E.Qt(day);
% Intermediary circulation
Qi = E.Qi(day);
% Length scale for velocity gradient
%DiffZ = 0.5*(D.H(2) + D.H(3));
if D.DWR == 0
    DiffZ = D.H(2);
else
    DiffZ = SillData.Hmax(1);
end
% Square of dz for velocity shear
dzsq=DiffZ^2;
% Density difference between layers
DelSig = D.rho(3) - D.rho(2); 
% Rough mean density in kgm-3
MeanSig = 0.5 * (D.rho(2) + D.rho(3)); 
% Approximate buoyancy frequency squared (s-2)
N_sq = (Const.g/MeanSig) * (DelSig/DiffZ); 
N_sq = max([N_sq 1e-6]);
% Length of loch
L = LochData.Len; 
% Get the x-sectional area of layer two
xarea_2 = D.xarea(2);
% Get the interfacial area between layers two and three
A_ifc = D.A_ifc(2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive estuarine induced flow in intermediate layer
U2e = (Qe + Qi) / xarea_2; 
% Amplitude of tidal velocity over the sill
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
    %disp(num2str([Param.t/86400 U2e U0t N_sq max(dudz2) min(Ri) mean(w)]))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entrainment velocity from layer 3 into layer 2
w23 = mean(w);

% if deep water renewal is occurring, entrainment is downwards
if D.DWR > 0; 
    % Scale the entrainment to excessive entrainment into layer 3 when
    % layer 3 is thin.
    Q_scale = (Const.H3max - D.H(3)) / Const.H3max;
    Q_scale = max([min([Q_scale 1]) 0]);
    w23 = -w23 * Q_scale; 
else
    Q_scale = D.H(3) / Const.H3max;
    w23 = w23 * Q_scale; 
end;

% convert to entrainment flux in m3s-1
% N.B. Vertical fluxes positive upwards
Ent23 = A_ifc * w23;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end