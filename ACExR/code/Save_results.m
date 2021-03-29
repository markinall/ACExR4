function save_results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function save_results
%
% Saves the final set of results to file. This routine truncates any
% spin-up period at the start of the output arrays.
%
%
% Phil Gillibrand
% 15/6/2007
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('Saving results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r1 and rn specify the first and last valid records in the output 
% arrays.
r1 = Param.ndays_spin_up + 1;
rn = Param.Ndays;

% Volume fluxes
E.Qe = E.Qe(r1:rn);
E.Qt = E.Qt(r1:rn);
E.Qi = E.Qi(r1:rn);
E.Qh = E.Qh(r1:rn,:);
E.Qw12 = E.Qw12(r1:rn);
E.Qw23 = E.Qw23(r1:rn);
E.Kz12 = E.Kz12(r1:rn);
E.Kz23 = E.Kz23(r1:rn);
E.Kzb = E.Kzb(r1:rn,:);
E.Qf = E.Qf(r1:rn);
E.Ue = E.Ue(r1:rn,:);
E.Ut = E.Ut(r1:rn,:);

% Layer parameters
Param.H = Param.H(r1:rn,:);
Param.S = Param.S(r1:rn,:);
Param.T = Param.T(r1:rn,:);
Param.V = Param.V(r1:rn,:);
Param.rho = Param.rho(r1:rn,:);
if Param.Ntracer > 0
	Param.C1 = Param.C1(r1:rn,:);
end
if Param.Ntracer > 1
	Param.C2 = Param.C2(r1:rn,:);
end
Param.N = Param.N(r1:rn,:);
Param.Ndays = size(Param.H,1);
Param.Kz = Param.Kz(r1:rn,:);
Param.Kzb = Param.Kzb(r1:rn,:);
Param.DWR = Param.DWR(r1:rn,:);
Param.Ri = Param.Ri(r1:rn,:);

% Boundary forcing
Bdata.Uw = Bdata.Uw(r1:rn);
Bdata.Qf = Bdata.Qf(r1:rn);
if Param.OBdata == 1
    Bdata.S_ext = Bdata.S_ext(:,r1:rn);
    Bdata.T_ext = Bdata.T_ext(:,r1:rn);
    Bdata.rho_ext = Bdata.rho_ext(:,r1:rn);
    if Param.Ntracer > 0
    	Bdata.C1_ext = Bdata.C1_ext(:,r1:rn);
    end
    if Param.Ntracer > 1
    	Bdata.C2_ext = Bdata.C2_ext(:,r1:rn);
    end
end
Bdata.a0 = Bdata.a0(r1:rn);
Bdata.Qshf = Bdata.Qshf(r1:rn);
Bdata.T0 = Bdata.T0(r1:rn,:);
Bdata.S0 = Bdata.S0(r1:rn,:);
Bdata.rho0 = Bdata.rho0(r1:rn,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save final set of results to file
save results E Param Bdata SillData LochData Const Hypso

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
