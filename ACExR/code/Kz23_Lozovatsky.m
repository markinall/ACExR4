function KZ = Kz23

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function KZ = Kz23
%
% evaluates definite symbolic integral of the Kz(x) along the axis
% of a basin. Kz is based on the Ri parameterisation of 
% Lozovatsky et al (DSR I, 2006), the x dependency comes through
% the tendency towards zero of the estuarine circulation as the 
% head of the loch is approached.
%
% called from CalcE function.
%
% for reference:
%
% A = kappa * ustar * H2 , where kappa is Von Karmans const 
% (0.4), ustar is the friction velocity, and H2 is the 
% intermediate layer thickness.
%
% B = 10 * Nsquared
%
% C = sqrt (U2e + U2t) , where U2e is the intermediate layer 
% estuarine velocity near the sill and U2t the tidal velocity; 
% the bottom layer is assumed to be motionless.
%
% Mark Inall
% 24 March 2006
% Modified by Paul Tett, November 2006
% Modified by Phil Gillibrand, June 2007.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get required parameters
day = Param.day;
% Estuarine exchange
Qe = abs(D.Qe); 
% Tidal exchange
Qt = E.Qt(day);
% Salinity difference between layers 2 and 3
DelS = D.S(3) - D.S(2); 

% Calculate KZ if the water column is stable, otherwise mixing is enhanced
% by overturning (water column instability).
if DelS > 0.0 
    % Distance between layer centres
    DiffZ = 0.5 * (D.H(2) + D.H(3));
    % Linear EoS conversion to density
    DelSig = Const.beta * DelS; 
    % Rough mean density in kgm-3
    MeanSig = 1025; 
    % Approximate buoyancy frequency squared (s-2)
    N_sq = (Const.g/MeanSig) * (DelSig/DiffZ); 
    % N^2 values constrained between +1e-6 and +1e-2 s-2
    N_sq = min([max([N_sq 1e-6]) 1e-2]);
    % Length of loch
    L = LochData.Len; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the x-sectional area of layer two
    xarea_2 = D.xarea(2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Derive estuarine induced flow in intermediate layer +ve
    U2e = max([Qe / xarea_2 0.01]); 
    % Derive tidally induced flow in intermediate layer +ve
    U2t = max([Qt / xarea_2 0.01]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Some constants for the Kz parameterisation - these are assumed 
    % and have been used to derive the Kz expression
    % RiCr = 0.1; % critical Richardson Number
    % P = 1; % tuning exponent
    % PrTr = 1 + Ri/RiB, RiB = 0.1 - using to determine turbulent 
    %        Prandtl Number (PrTr)
    % Von Karmans constant
    kappa = Const.kappa; 
    % Friction velocity
    ustar = sqrt(Const.Cd) * U2t;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some constant combinations for a given time step
    % A = kappa * ustar * H(day,2); % surface layer formulation
    % Pacanowski and Philander constant background diffusivity
    A = 1e-4;
    B = 10 * N_sq * DiffZ * DiffZ;
    C = sqrt(U2t + U2e);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %expression for Kz from Lozovatsky et al (DSR I, 2006) - integrated along axis
    KZ = -1/2*A*(-4*(C^2+C*B)^(1/2)*C^2-6*(C^2+C*B)^(1/2)*C*B- ...
          3*B^2*(C^2+C*B)^(1/2)+10*B^2*atan(C/((C+B)*C)^(1/2))*C+ ...
          3*B^3*atan(C/((C+B)*C)^(1/2))+8*B*atan(C/((C+B)*C)^(1/2))*C^2)...
          /(C+B)/(2*C+B)/(C^2+C*B)^(1/2);
    %PT == KZ in m2/s ??
    %PT == minimum will always be 1e-4 
    %PT == (1e-5 option not called because N_sq made always +ve)
    %PT == implies minimum 10 m2/d and timescale of 4h(3)^2/KZ
    %PT == about 40 days for h(3) = 10 m -- Creran
    %PT == about 1000 days for h(3) = 50 m -- Etive
    %PT == look at Edwards, Grantham estimates of Etive Kz
    %PT == ----------------------------------------------------------
else 
    %PT == end of Mark's calculation of Kz, for DelS +ve
    %PT == alternative calculation for overturn, with constraints
    MaxKz=D.H(3)^2/(2*86400); %PT == constraint on maximum mixing velocity
	%PT == about 5e-4 for Creran
    MinKz=D.H(3)^2/(7*86400);%PT == ensures overturn within week
    AltKz=1e-2;%PT == proposed default value for overturn
    KZ=min(MaxKz, min(AltKz, MinKz));%PT == m2/s
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate interfacial areas between layers
A_23 = D.A_ifc(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert to a flux of m3/s. 
H23 = 0.5 * (D.H(2) + D.H(3));
KZ = KZ * A_23 / H23;

% end function
end
