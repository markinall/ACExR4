function eps = Efficiency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function eps = Efficiency
%
% Estimates the efficiency of tidal exchange. Based on the work of
% Wells & van Heijst (2003) who extended the work of Stommel and
% Farmer (1952). The efficiency calculation is based on the jet-sink 
% model of estuarine exchange including the dipole effect of 
% Wells & van Heijst (2003, Dyn. Atmos. Ocean., 37, 223-244).
%
% Usage:    LochData contains loch information inc. tidal range
%           SillData contains topographic data for each sill
%           Param contains various parameters including layer depths
%           Const contains tidal period
%
% Phil Gillibrand
% 27/3/2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define topographic parameter in calculation
if LochData.Nsill == 0
    U = LochData.Range * Param.Af / Const.Tperiod;
elseif SillData.Current(1) > 0
    %U = 0.5 * SillData.Current(1);
    U = SillData.Current(1);
else
    U = 2 * (LochData.HWarea + LochData.LWarea) * LochData.Range / ...
        (SillData.Xarea(1) * Const.Tperiod);
end
WoUT = Param.Bm / (U * Const.Tperiod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate exchange efficiency
eps = 1 - sqrt(4 * WoUT / pi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include dipole escape factor
if WoUT < 0.13
    eps = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default value as backup for missing values
if ~isfinite(eps)
    eps = 0.5;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Efficiency should be applied on both the flood and the ebb tides.
eps = eps^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end