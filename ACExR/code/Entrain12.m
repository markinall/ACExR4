function Ent12 = Entrain12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ent12 = Entrain12
% numerical evaluation of entrainment from layer 1 into layer 2
% by barotropic tidal current subducted under layer 1
% Based on the Ri parameterisation of Ellison and Turner
%
% called from CalcE function.
%
% for reference: A, B, and C are tunable constants
%
A=0.00008;,B=0.00;,C=10;
% Mark Inall
% 29 March 2006

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Salinity difference between layers 1 and 2.
DelS = D.S(2) - D.S(1); 
% Linear EoS conversion to density
DelSig = Const.beta * DelS; 
% Rough mean density in kgm-3
MeanSig = 1025; 
% Approximate buoyancy frequency squared (s-2)
N_sq = (Const.g/MeanSig) * (DelSig/DiffZ); 
N_sq = min([max([N_sq 1e-6]) 1e-2]);
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
    U0t = Qt / SillData.Xarea(1);
else
    U0t = Qt / (Hypso.vol(end) / LochData.Len);
end
% Mean tidal velocity in layer 2
U2t = Qt / xarea_2;
% Length scale of velocity decay in layer 2
Xl = L;
if U0t > 0; Xl = L * U2t / U0t; end
% M2 Tidal period
T=Const.Tperiod;
omega = 2*pi/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate entrainment for the period from the present time step 
% to the next time step.
t = Param.t;
x=0:20:L-10;
w=zeros(1,length(x));
for j = 1:length(x)
    U2 = U2e + U0t*sin(omega*t)*exp(-x(j)/Xl);
    if U2 > 0
        dudz2 = (U2*U2/dzsq);
        Ri = max([N_sq/dudz2 0.1]);
        w(j) = (A-B*Ri)/(1+C*Ri);
    end
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
