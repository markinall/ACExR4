function [Ynew] = tracer(Y,Yf,Y0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Ynew] = tracer(Y,Yf,Y0)
%
% Calculates values of tracer in each layer (Y(1), Y(2), Y(3)) for 
% day number 'day' forced by volume exchanges in each layer.
%
% Usage:    Y contains the present values of the parameter to be 
%                      determined e.g. salinity
%           Yf is the river flow boundary value of Y
%           Y0 is the oceanic boundary value of Y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define local variables
% Day
day = Param.day;
% Layer volumes at previous time step
Vp = D.Vp;
% Layer volumes
V = D.V;
% Most recent flux values
Qt = E.Qt(day);
Qi = E.Qi(day);
Qf = E.Qf(day);
Qe = D.Qe;
Qh = D.Qh;
Qw12 = D.Qw12;
Qw23 = D.Qw23;
Kz12 = D.Kz12;
Kz23 = D.Kz23;

% Initialise output array
Ynew = zeros(1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 1
Ynew(1) = Y(1) * Vp(1) / V(1) ...
        + Const.deltaT * ((1/V(1)) ...
        * (Qe * (Y(2) - Y(1)) ...
        +  Qf * (Yf - Y(1)) ...
        +  Qh * Y(2) ...
        +  Qw12 * Y(1) ...
        +  Kz12 * (Y(2) - Y(1))) ...
        + (1/(V(1) + V(2))) * Qi * (Y0(1) - Y(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 2
Ynew(2) = Y(2) * Vp(2) / V(2) ...
        + Const.deltaT * ((1/V(2)) ...
        * ((1 - D.DWR) * Qe * (Y0(2) - Y(2)) ...
        -  D.DWR * Qe * Y(2) ...
        -  Qh * Y(2) ...
        -  Qw12 * Y(1) ...
        -  Kz12 * (Y(2) - Y(1)) ...
        +  E.Qt(day) * (Y0(2) - Y(2))) ...
        + (1/(V(1) + V(2))) * Qi * (Y0(2) - Y(2)) ...
        + Const.relaxT * (Y0(2) - Y(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 3 (and fluxes between layers 2 and 3)
if LochData.Nsill > 0
    Ynew(2) = Ynew(2) + Const.deltaT * (1/V(2)) ...
            * (0.5 * (Qw23 + abs(Qw23)) * Y(3) ...
            +  0.5 * (Qw23 - abs(Qw23)) * Y(2) ...
            +  Kz23 * (Y(3) - Y(2)));
    Ynew(3) = Y(3) * (Vp(3) / V(3)) ...
            + (Const.deltaT / V(3)) ...
            * (D.DWR * Qe * Y0(2) ...
            -  0.5 * (Qw23 + abs(Qw23)) * Y(3) ...
            -  0.5 * (Qw23 - abs(Qw23)) * Y(2) ...
            +  Kz23 * (Y(3) - Y(2)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
