function Hypsography

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Hypsography
%
% Derives hypsographic function for sea loch
% 
% Usage:    LochData contains data from the sea loch catalogue
%           Hypso contains a profile of horizontal planar areas
%               (A), a cumulative vertical cross-sectional area
%               (xarea), cumulative layer volume (vol) and 
%               the depths of the planar areas (z).
%           
% This version based on sea loch catalogue data vlaues assuming 
% an inverted pyramid shape. All variables start at z = 0, so 
% the first values of vol and xarea are zero.
%
% Phil Gillibrand
% 22/3/2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('Deriving hypsographic function');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get basic topographic data from catalogue-derived array. If
% the "entrance" sill is designated as an internal sill, first
% adjust all the parameters in LochData and SillData accordingly.
% The adjustment is based on the ratio between the LW surface area
% upstream of the sill and the whole loch LW surface area.
if Param.Esill > 1
    isill = Param.Esill;
    afac = SillData.LWuparea(isill) / LochData.LWarea;
    LochData.Len = LochData.Len * afac;
    LochData.HWarea = SillData.HWuparea(isill);
    LochData.LWarea = SillData.LWuparea(isill);
    LochData.Area2m = round(LochData.Area2m * afac);
    LochData.Area5m = round(LochData.Area5m * afac);
    LochData.Area10m = round(LochData.Area10m * afac);    
    LochData.LWvol = round(LochData.LWvol * afac);
    LochData.HWvol = round(LochData.HWvol * afac);
    LochData.Wshed = round(LochData.Wshed * afac);
    LochData.Qf = round(LochData.Qf * afac);
    LochData.Hmean = round(LochData.LWvol *10 / LochData.LWarea) / 10 ;
    LochData.Nsill = LochData.Nsill - isill + 1;
    % Correct the mean depth to MSL
    LochData.Hmean = LochData.Hmean + LochData.Range / 2;
    % Re-calculate flushing time
    LochData.Tf = 89424 * LochData.LWvol / ...
        ((LochData.HWarea + LochData.LWarea)*0.7*LochData.Range);
    % Remove the data for sills outside the new "entrance" sill
    SillData.Len = SillData.Len(isill:end,:);
    SillData.HWwid = SillData.HWwid(isill:end,:);
    SillData.LWwid = SillData.LWwid(isill:end,:);
    SillData.Hmax = SillData.Hmax(isill:end,:);
    SillData.Hmean = SillData.Hmean(isill:end,:);
    SillData.Xarea = SillData.Xarea(isill:end,:);
    SillData.HWuparea = SillData.HWuparea(isill:end,:);
    SillData.LWuparea = SillData.LWuparea(isill:end,:);
    SillData.Current = SillData.Current(isill:end,:);
    SillData.Hbasin = SillData.Hbasin(isill:end,:);
    SillData.Barea = SillData.Barea(isill:end,:);
    Param.Bm = SillData.Xarea(1) / SillData.Hmean(1);
end

% Set key parameters for hysography
% Surface area at mean sea level.
Af = 0.5 * (LochData.HWarea + LochData.LWarea);
% Estuary length
L = LochData.Len;

% The maximum depth is set so that the total volume and surface area of the
% estuary match those of the real system when the estuary is assumed to
% have the shape of an inverted pyramid.
% Calculate the estuary volume below mean sea level
Vol = 0.5 * (LochData.HWvol + LochData.LWvol);
LochData.Hmax = 3 * Vol / Af;
Hmax = LochData.Hmax;
SillData.Hbasin = min(SillData.Hbasin, Hmax);
if LochData.Hmax < 0; return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate hysographic function
Hypso.dz = 1; 
if LochData.Hmax < 10; Hypso.dz = 0.1; end
dzr = round(1/Hypso.dz);
LochData.Hmax = floor(LochData.Hmax / Hypso.dz) / dzr;
Hmax = LochData.Hmax;
Hypso.z = 0:Hypso.dz:Hmax;
Hypso.A = (Hmax - Hypso.z).^2 * Af / (Hmax^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate cumulative cross-sectional area
% First get layer volumes
Nlayer = length(Hypso.z);
for ic = 2:Nlayer
    LayerVol(ic) = (Hypso.A(ic) + Hypso.A(ic-1)) ...
                 * (Hypso.z(ic) - Hypso.z(ic-1)) ...
                 / 2;
end
%Xarea = LayerVol * Hmax ./ ((Hmax - Hypso.z) * L);
Xarea = LayerVol / L;
% Specify the cross-sectional area as a cumulative sum
Hypso.xarea = cumsum(Xarea);
% Store the cumulative volume
Hypso.vol = cumsum(LayerVol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

end
