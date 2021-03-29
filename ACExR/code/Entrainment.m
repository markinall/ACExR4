function Ent = Entrainment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ent = Entrainment
% numerical evaluation of entrainment from layer 1 into layer 2
% by barotropic tidal current subducted under layer 1
% Based on the Ri parameterisation of Ellison and Turner
%
% called from CalcE function.
%
% for reference: A, B, and C are tunable constants for the Turner scheme.
%
% A=0.08;,B=0.1;C=5;
% Mark Inall
% 29 March 2006

%%%%%%%%%% MEI Comment 14th Dec 2007 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error corrected in orogonal code - factor of U missing in entrainment
% velocity. the code now runs with the A,B,C values as given by Turner
%correct values from Turner 1986 are: A=0.08;,B=0.1;,C=5; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEI Comment - Feb 2008 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entrainment param of Princevac et al, JFM 2005, vol533 pp259-268
% implemented - Turner 1986 param no longer used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MEI - 12 Feb 2008 comment
% shear between L1 and L2 now includes the estuarine velocity of the upper
% layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAG - 12 July 2011
% Calculates entrainment across BOTH horizontal interfaces (i.e. combines 
% former Entrain12_inall and Entrain23_inall scripts into one)
% using either Turner, Princevac et al, or Strang & Fernando, JFM 2001, 
% vol 428, pp349-386, algorithms. Entrainment is calculated using a sub-
% time-stepping algorithm to improve numerical stability.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ent = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% No entrainment is modelled if deep water renewal is occurring
%if D.DWR == 1
    %return
%end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Otherwise, model entrainment.
% Get required parameters
day = Param.day;
% Temporary layer thicknesses
H_tmp = D.H;
% Estuarine exchange
Qe = D.Qe; 
% Tidal exchange
Qt = E.Qt(day);
% Intermediary circulation
Qi = E.Qi(day);

%% Derive stratification parameters for Layer 1 - Layer 2 interface
% Vertical distance between layer centres
DiffZ_1 = 0.5*(D.H(1) + D.H(2));
DiffZ_1 = D.H(2);
% Square of dz for velocity shear
dzsq=D.H(2)^2;
% Density difference between layers 1 and 2.
DelSig_1 = D.rho(2) - D.rho(1);
% Rough mean density in kgm-3
MeanSig_1 = 0.5 * (D.rho(1) + D.rho(2)); 
% Buoyancy difference
DelB_1 = Const.g * DelSig_1 / MeanSig_1;
% Approximate buoyancy frequency squared (s-2)
N_sq_1 = DelB_1 / DiffZ_1; 
N_sq_1 = max([N_sq_1 1e-6]);

%% Derive stratification parameters for Layer 2 - Layer 3 interface
if LochData.Nsill > 0
    DiffZ_2 = D.H(2);
    % Square of dz for velocity shear
    % dzsq_2=DiffZ_2^2;
    % Density difference between layers
    DelSig_2 = D.rho(3) - D.rho(2); 
    % Rough mean density in kgm-3
    MeanSig_2 = 0.5 * (D.rho(2) + D.rho(3)); 
    % Buoyancy difference
    DelB_2 = Const.g * DelSig_2 / MeanSig_2;
    % Approximate buoyancy frequency squared (s-2)
    N_sq_2 = DelB_2 / DiffZ_2; 
    N_sq_2 = max([N_sq_2 1e-6]);
end

%% Specify remaining parameters
% Length of loch
L = LochData.Len; 
% Get the x-sectional area of layer one
xarea_1 = D.xarea(1);
% Get the x-sectional area of layer two
xarea_2 = D.xarea(2);

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
%Length scale of velocity decay
Xl = L/5;
% M2 Tidal period
T=Const.Tperiod;
omega = 2*pi/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate entrainment across interface over the duration
%% of the current time step.
w12 = 0; w23 = 0;
Ri = 0;
x=0:20:L-10;
dt = 120;
N_substeps = Const.deltaT / dt;
t = Param.t - Const.deltaT;
for ic = 1:N_substeps
    t = t + dt;
    w12_tmp=zeros(1,length(x));
    w23_tmp=zeros(1,length(x));
    U2 = U2e + U2t*sin(omega*t);
    if U2 > 0
        U2x = U2*exp(-x/Xl);
        %% Entrainment across upper interface
        if D.DWR == 0
            dudz2 = (U2x .* U2x / dzsq);
            Ri = N_sq_1 ./ dudz2;
    
            % Turner (1986) formulation
            % A=0.08;,B=0.1;C=5;
            %id = Ri < 0.8;
            %w12_tmp(id) = U2x(id) .* (A - B * Ri(id)) ./ (1 + C * Ri(id));
        
            % Princevac et al (2005) formulation
             Ri = max(Ri, 0.15);
             w12_tmp = U2x .* (0.05 * Ri.^(-0.75));
        
            % Strang & Fernando (2001) formulation, using Region I and III eqns
            % Ri = DelB_1 * H_tmp(2) ./ (U2x .* U2x);
            % Ri(Ri < 0.87) = 0.87;
            % w12_tmp = 0.02 * U2x .* Ri .^ (-1.3);
        
            % Scale the entrainment to prevent loss of surface layer volume
            % when surface layer is already at its minimum value.
            Q_scale = (4*(((Const.H1max - H_tmp(1)) * (H_tmp(1) - Const.H1min)) / ...
                ((Const.H1max - Const.H1min) ^ 2))).^2;
            Q_scale = max([min([Q_scale 1]) 0]);
            % N.B. Vertical fluxes positive upwards
            w12_tmp = -w12_tmp * Q_scale;
        end
        
        %% Entrainment across lower interface
        if LochData.Nsill > 0
            dudz2 = (U2x .* U2x / dzsq);
            Ri = N_sq_2 ./ dudz2;
            
            % Turner (1986) formulation
            % A=0.08;,B=0.1;C=5;
            % id = Ri < 0.8;
            % w23_tmp(id) = U2x(id) .* (A - B * Ri(id)) ./ (1 + C * Ri(id));
            
            % Princevac et al (2005) formulation
            Ri = max(Ri, 0.15);
            w23_tmp = U2x .* (0.05 * Ri.^(-0.75));
        
            % Strang & Fernando (2001) formulation, using Region I and III eqns
            % Ri = DelB_2 * H_tmp(2) ./ (U2x .* U2x);
            % Ri(Ri < 0.87) = 0.87;
            % w23_tmp = 0.02 * U2x .* Ri .^ (-1.3);
            
            % Scale the entrainment to prevent loss of surface layer volume
            % when bottom layer is already at its minimum value.
            if D.DWR == 0
                Q_scale = ((H_tmp(3) - Const.H3min) / (Const.H3max - Const.H3min)).^2;
            else
                Q_scale = (4*(((Const.H3max - H_tmp(3)) * (H_tmp(3) - Const.H3min)) / ...
                    ((Const.H3max - Const.H3min) ^ 2))).^2;
                %Q_scale = ((Const.H3max - H_tmp(3)) / (Const.H3max - Const.H3min)).^2;
                w23_tmp = -w23_tmp;
            end
            Q_scale = max([min([Q_scale 1]) 0]);
            w23_tmp = w23_tmp * Q_scale;
        end
    end
    
    % Get mean entrainment along length of loch basin
    w12_tmp = mean(w12_tmp);
    w23_tmp = mean(w23_tmp);
    
    % Recalculate temporary layer thickness & cross-section
    H_new(1) = H_tmp(1) + dt * w12_tmp;
    H_new(2) = H_tmp(2) + dt * (w23_tmp - w12_tmp);
    H_new(3) = H_tmp(3) - dt * w23_tmp;
    if H_new(1) < Const.H1min
        w12_tmp = (Const.H1min - H_tmp(1)) / dt;
        H_new(1) = H_tmp(1) + dt * w12_tmp;
        H_new(2) = H_tmp(2) + dt * (w23_tmp - w12_tmp);
    end
    if H_new(1) > Const.H1max
        w12_tmp = (Const.H1max - H_tmp(1)) / dt;
        H_new(1) = H_tmp(1) + dt * w12_tmp;
        H_new(2) = H_tmp(2) + dt * (w23_tmp - w12_tmp);
    end
    %disp(num2str(H_tmp))
    
    % Update temporary layer thicknesses
    H_tmp = H_new;
    
    % Sum entrainment rates to get temporal mean later
    w12 = w12 + w12_tmp;
    w23 = w23 + w23_tmp;
    %disp(num2str([H_tmp Qe Qt max(U2x) DelB_1 min(Ri) mean(w12)]))

    % Recalculate cross-sectional areas
    h2 = min([H_tmp(1)+H_tmp(2), Hypso.z(end)]);
    xarea_1 = getX(H_tmp(1));
    xarea_2 = getX(h2) - xarea_1;
    
    % Recalculate shear length scale
    dzsq = H_tmp(2)^2;

    % Recalculate buoyancy frequency
    N_sq_1 = DelB_1 / H_tmp(2); 
    N_sq_1 = max([N_sq_1 1e-6]);
    if LochData.Nsill > 0
        N_sq_2 = DelB_2 / H_tmp(2); 
        N_sq_2 = max([N_sq_2 1e-6]);
    end

    % Recalculate layer velocities
    U2t = pi * Qt / (2 * Param.eps * xarea_2);
    U2e = (Qe + Qi) / xarea_2; 
end
    
% Get time mean for w12
w12 = w12 / N_substeps;
w23 = w23 / N_substeps;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provide entrainment flux in m3s-1
%disp(num2str([D.H D.DWR Q_scale w12 0.5*Const.deltaT*w12]))
A_ifc = getA(D.H(1) + 0.5 * Const.deltaT * w12);
Ent(1) = A_ifc * w12;
if LochData.Nsill > 0
    A_ifc = getA(D.H(1) + D.H(2) + 0.5 * Const.deltaT * w23);
    Ent(2) = A_ifc * w23;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = getX(H)

global Hypso

iz1 = floor(H / Hypso.dz) + 1;
iz2 = ceil(H / Hypso.dz) + 1;
f = (H - Hypso.z(iz1)) / Hypso.dz;
X = (1-f)*Hypso.xarea(iz1) + f*Hypso.xarea(iz2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
