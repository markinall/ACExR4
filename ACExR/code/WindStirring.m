function WE = WindStirring
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function WE = WindStirring
% 
% Calculation of wind-driven entrainment (stirring), based on Stigenrandt
% (1985).
% called from CalcE function.
% Phil Gillibrand
% June 2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%Initialise WE
WE = [0 0];

% set day counter
day = Param.day;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate entrainment velocity,w_e, from Layer 2 to Layer 1 due to wind
% stirring. See eqns (6), (8) and (9) from Stigebrandt (1985).
% Entrainment velocity is positive upwards.
% First, calculate (and limit) reduced gravity
drho = max([D.rho(2) - D.rho(1), 0.1]);
g_red = Const.g * drho / D.rho(2);

% Calculate the total buoyancy input B (m2/s3) to the surface layer
% from heat and freshwater inputs (Stigebrandt, 1985).
B = (Const.g * Const.alpha * D.Qshf / (D.rho(1)^2 * Const.cp)) ...
    + (Param.M * E.Qf(day) / Param.Af);
if B < 0
    epsilon = 0.05;
else
    epsilon = 1;
end
        
% Calculate the entrainment velocity (m/s, +ve vertically upwards)
w_e = (Param.N(day) / (g_red * D.H(1))) - (epsilon * B / g_red);
     
if w_e > 0 
    % Calculate entrainment velocity with fourth-order Runge-Kutta
    h1 = D.H(1) + Const.deltaT * w_e;
    w_e = (Param.N(day) / (g_red * h1)) - (epsilon * B / g_red);
    h2 = D.H(1) + 0.5 * Const.deltaT * w_e;
    w_e = (Param.N(day) / (g_red * h2)) - (epsilon * B / g_red);
    h3 = D.H(1) + 0.5 * Const.deltaT * w_e;
    w_e = (Param.N(day) / (g_red * h3)) - (epsilon * B / g_red);
    h4 = D.H(1) + Const.deltaT * w_e;
    dh = (h1 + 2*h2 + 2*h3 + h4)/6;
    % Constrain updated surface layer thickness
    if dh > Const.H1max; dh = Const.H1max; end
    % derive entrainment velocity
    w_e = (dh - D.H(1)) / Const.deltaT;
end
   
if w_e < 0
    % In a retreating regime (i.e. w_e < 0), the pycnocline retreats to a 
    % depth given by the Monin-Obukhov length (Stigebrandt, 1985): 
    % h = Param.N / B
    % We assume that this occurs over an e-folding time scale of 24
    % hours, giving:
    H_MO = Param.N(day) / B;
    % Constrain by H1 limits
    H_MO = min(Const.H1max, max(H_MO, Const.H1min));
    w_e = (H_MO - D.H(1)) / Const.Tday;
end

% Specify the volume flux due to entrainment
A_ifc = getA(D.H(1) + 0.5 * Const.deltaT * w_e);

% Calculate the volume flux
WE(1) = A_ifc * w_e;

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if LochData.Nsill > 0
    %% Calculate entrainment velocity,w_e, from Layer 3 to Layer 2 due to wind
    % stirring. See eqns (6), (8) and (9) from Stigebrandt (1985).
    % Entrainment velocity is positive upwards.
    % First, calculate (and limit) reduced gravity
    drho = max([D.rho(3) - D.rho(2), 0.1]);
    g_red = Const.g * drho / D.rho(3);
    
    % Calculate the entrainment velocity (m/s, +ve vertically upwards)
    h0 = D.H(1) + D.H(2);
    w_e = (Param.N(day) / (g_red * h0)) - (epsilon * B / g_red);
    
    if w_e > 0 
        % Calculate entrainment velocity with fourth-order Runge-Kutta
        h1 = h0 + Const.deltaT * w_e;
        w_e = (Param.N(day) / (g_red * h1)) - (epsilon * B / g_red);
        h2 = h0 + 0.5 * Const.deltaT * w_e;
        w_e = (Param.N(day) / (g_red * h2)) - (epsilon * B / g_red);
        h3 = h0 + 0.5 * Const.deltaT * w_e;
        w_e = (Param.N(day) / (g_red * h3)) - (epsilon * B / g_red);
        h4 = h0 + Const.deltaT * w_e;
        dh = (h1 + 2*h2 + 2*h3 + h4)/6;
        % Constrain updated surface layer thickness
        dh_max = LochData.Hmax - Const.H3min;
        if dh > dh_max; dh = dh_max; end
        % derive entrainment velocity
        w_e = (dh - h0) / Const.deltaT;
    end
    
    if w_e < 0
        % In a retreating regime (i.e. w_e < 0), the pycnocline retreats to a 
        % depth given by the Monin-Obukhov length (Stigebrandt, 1985): 
        % h = Param.N / B
        % We assume that this occurs over an e-folding time scale of 24
        % hours, giving:
        H_MO = Param.N(day) / B;
        % Constrain by H3 limits
        h0min = LochData.Hmax - Const.H3max;
        H_MO = max(H_MO, h0min);
        w_e = (H_MO - h0) / Const.Tday;
    end
    
    % Specify the volume flux due to entrainment
    A_ifc = getA(h0 + 0.5 * Const.deltaT * w_e);
    
    % Calculate the volume flux
    WE(2) = A_ifc * w_e;
    WE(2)=0;
end
    
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
