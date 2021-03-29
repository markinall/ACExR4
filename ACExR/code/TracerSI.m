function [Ynew] = tracerSI(Y,Yf,Y0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Ynew] = tracerSI(Y,Yf,Y0)
%
% Calculates values of tracer in each layer (Y(1), Y(2), Y(3)) for 
% day number 'day' forced by volume exchanges in each layer.
%
% Usage:    Y contains the present values of the parameter to be 
%                      determined e.g. salinity
%           Yf is the river flow boundary value of Y
%           Y0 is the oceanic boundary value of Y
%
% Note:     This version solves the vertical diffusion terms implicitly,
%            while the remaining terms are still solved explicitly.
%
% Philip Gillibrand
% 18 February 2009.
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

% Initialise output array
Ynew = zeros(1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 1
Ynew(1) = Y(1) * Vp(1) / V(1) ...
        + Const.deltaT * ((1/V(1)) ...
        * (Qe * (Y(2) - Y(1)) ...
        +  Qf * (Yf - Y(1)) ...
        +  0.5 * (Qh(1) + abs(Qh(1))) * Y(2) ...
        +  0.5 * (Qh(1) - abs(Qh(1))) * Y(1) ...
        +  Qw12 * Y(1)) ...
        + (1/(V(1) + V(2))) * Qi * (Y0(1) - Y(1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 2
Ynew(2) = Y(2) * Vp(2) / V(2) ...
        + Const.deltaT * ((1/V(2)) ...
        * ((1 - D.DWR) * (Qe + Qt) * (Y0(2) - Y(2)) ...
        +  D.DWR * Qe * (Y(3) - Y(2)) ...
        -  0.5 * (Qh(1) + abs(Qh(1))) * Y(2) ...
        -  0.5 * (Qh(1) - abs(Qh(1))) * Y(1) ...
        +  0.5 * (Qh(2) + abs(Qh(2))) * Y(3) ...
        +  0.5 * (Qh(2) - abs(Qh(2))) * Y(2) ...
        -  Qw12 * Y(1)) ...
        + (1/(V(1) + V(2))) * Qi * (Y0(2) - Y(2)) ...
        + Const.relaxT * (Y0(2) - Y(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Layer 3 (and fluxes between layers 2 and 3)
if LochData.Nsill > 0
    Ynew(2) = Ynew(2) + Const.deltaT * (1/V(2)) ...
            * (0.5 * (Qw23 + abs(Qw23)) * Y(3) ...
            +  0.5 * (Qw23 - abs(Qw23)) * Y(2));
    Ynew(3) = Y(3) * (Vp(3) / V(3)) ...
            + (Const.deltaT / V(3)) ...
            * (D.DWR * (Qe + Qt) * (Y0(2) - Y(3)) ...
            -  0.5 * (Qh(2) + abs(Qh(2))) * Y(3) ...
            -  0.5 * (Qh(2) - abs(Qh(2))) * Y(2) ...
            -  0.5 * (Qw23 + abs(Qw23)) * Y(3) ...
            -  0.5 * (Qw23 - abs(Qw23)) * Y(2));
    %Call function for implicit diffusion
    Ynew = implicit_diffusion(Ynew);
else
    Ynew(1:2) = implicit_diffusion(Ynew(1:2));    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Yupdated = implicit_diffusion(Y)
%
% Solves vertical diffusion implicitly i.e.
% AY(n+1) = BY(n)
% (where Y(n) indicates the value of vector Y at time step n).
% so
% Y(n+1) = A'BY(n)
% where A' is the inverse of matrix A.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global Const D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = Const.deltaT;
Kz12 = D.Kz12;
Kz23 = D.Kz23;
V = D.V;
A = [];
B = [];

if length(Y) == 2
    % Build up matrices A and B.
    ilayer = 1;
    cjp1 = (dt * Kz12) / (2 * V(ilayer));
    A(ilayer,1:2) = [1+cjp1 -cjp1];
    B(ilayer,1:2) = [1-cjp1 cjp1];    
    ilayer = 2;
    cj = (dt * Kz12) / (2 * V(ilayer));
    A(ilayer,1:2) = [-cj 1+cj];
    B(ilayer,1:2) = [cj 1-cj];
else
    % Build up matrices A and B.
    ilayer = 1;
    cjp1 = (dt * Kz12) / (2 * V(ilayer));
    A(ilayer,1:3) = [1+cjp1 -cjp1 0];
    B(ilayer,1:3) = [1-cjp1 cjp1 0];    
    ilayer = 2;
    cj = (dt * Kz12) / (2 * V(ilayer));
    cjp1 = (dt * Kz23) / (2 * V(ilayer));
    A(ilayer,1:3) = [-cj 1+cj+cjp1 -cjp1];
    B(ilayer,1:3) = [cj 1-cj-cjp1 cjp1];
    ilayer = 3;
    cj = (dt * Kz23) / (2 * V(ilayer));
    A(ilayer,1:3) = [0 -cj 1+cj];
    B(ilayer,1:3) = [0 cj 1-cj];
end
     
    % get inverse of A
    Ainv = inv(A);

    % Implicit solution
    Yupdated = [Ainv * B * Y']';

end
