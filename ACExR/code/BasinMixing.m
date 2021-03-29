function W = BasinMixing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function W = BasinMixing
%
% Parameterisation of basin water mxing in Scottish
% sea lochs for deoxygenation risk assessment.
%
% Usage:    LochData contains loch information from catalogue
%           SillData contains data on each sill
%           Param contains mixed layer thickness and salinity
%           Const contains constants
%           R returns the rate of density reduction of basin
%               water (kg/m3/s)
%           W returns the work done in mixing the basin water
%           Te is the time scale (days) for renewal
%
% Phil Gillibrand
% July 2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Document Modifications Below
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp(' ');
%disp('Calculating basin mixing');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

% Initialise output variables
W = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Specify constants and other parameters
% Tidal (M2) frequency
Tm2 = Const.Tperiod;
omega = 2 * pi / Tm2;
% Reference density
rho0 = 1025;
% Flux Richardson Number (default)
Rf = 0.05;
% FjordEnv constant
Cw = 2;
% Number of sills in loch
Nsill = LochData.Nsill;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get values of mixed layer parameters
day = max([Param.day - 1, 1]);
% Density
rho = Param.rho(day,:);
% Mixed layer thickness
Hlay = Param.H(day,:);
% Tidal amplitude
a0 = Bdata.a0(Param.day);
% Surface area of loch basins
Basinarea = SillData.Barea;
% Upstream area at each sill
Aupstream = 0.5 * (SillData.HWuparea + SillData.LWuparea);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate time scale for renewal events
% Reduction in density required for deep-water renewal to occur
% Initial estimates from published literature (e.g. Edwards & 
% Edelsten, 1977; Gillibrand et al., 1995, 1996; Allen & 
% Simpson, 1998) suggest Re ~ 0.5 kg/m3
% The density reduction required in the first basin is determined
% by the standard deviation of the external density.
Ht = round(SillData.Hmax(1));
Re = std(Bdata.rho_ext(Ht,:));
% The following adjustment of Re for internal sills is ENTIRELY 
% ARBITRARY AND MUST BE IMPROVED.
for isill = 1:Nsill
    Re(isill) = Re(1) * log10(isill + 10);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Internal wave energy works against buoyancy in the deep layer.
% Waves are generated at the crest of sills in the stratified 
% tidal flow. Calculations of internal wave phase speed and energy
% flux are taken from Stigebrandt (1976, 1999, 2001) etc.
% Calculations are made for each sill, with waves assumed to be
% generated in both directions. The total energy flux is calculated
% and assumed to work against buoyancy in the water column below the
% shallowest sill depth. 
% The efficiency of mixing, the flux Richardson Number Rf, is 
% determined according to the Froude Number of the sill flow. If
% the flow is sub-critical, Rf = 0.05, if supercritical Rf = 0.01
% (Stigebrandt & Aure, 1989).

for isill = 1:Nsill
    
    % Work done against buoyancy in each basin
    W_basin = 0;
    
    % Set layer depths Ht and Hb, depending on sill and basin depths.
    Ht1 = max([SillData.Hmax(isill) Const.H1min]);
    Hb = floor(SillData.Hbasin(isill));
    
    % Following lines commented out on 8 November 2011.
    % If the surface layer is greater than sill depth, set
    % Ht = Hlay(1)
    %if Hlay(1) > Ht1
    %    Ht1 = Hlay(1);
    %end
    
    % Following lines inserted on 8 November 2011.
    if Hlay(1) < Ht1
        Ht1 = Hlay(1);
    end
    
    % Get the mean depth of the basin water
    uplayer = floor(Ht1) + 1;
    dnlayer = ceil(Ht1) + 1;
    f = Ht1 - uplayer + 1;
    At_loch = (1-f)*Hypso.A(uplayer) + f*Hypso.A(dnlayer);
    Vt_loch = (1-f)*Hypso.vol(uplayer) + f*Hypso.vol(dnlayer);
    uplayer = floor(Hb) + 1;
    dnlayer = ceil(Hb) + 1;
    f = Hb - uplayer + 1;
    Vb_loch = (1-f)*Hypso.vol(uplayer) + f*Hypso.vol(dnlayer) ...
                          - Vt_loch;
    if Basinarea(isill) > 0
        At = At_loch * Basinarea(isill) / Hypso.A(1);
        Vb = Vb_loch * Basinarea(isill) / Hypso.A(1);    
    else
        % otherwise use whole loch values to get mean depth
        At = At_loch;
        Vb = Vb_loch;
    end
    Hb1 = Vb / At;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Specify parameters needed for energy flux calculation
    % Reduced gravity
    delta_rho = rho(3) - rho(2);
    if Hlay(1) == Ht1
        delta_rho = rho(2) - rho(1);
    end
    %delta_rho = max([delta_rho 0.1]);
    if delta_rho < 0.1
        W(isill) = 0;
        continue
    end
    gprime = delta_rho * Const.g / rho(3);

    % calculate the internal wave phase speed, ci 
    ci = sqrt(gprime * Ht1 * Hb1 / (Ht1 + Hb1));
    % Calculate the energy flux due to internal tide propagating
    % upstream from the sill.
    Ej = 0.5 * rho0 * omega^2 * a0^2 * ci * (Hb1 / (Ht1 + Hb1)) ...
         * Aupstream(isill)^2 / SillData.Xarea(isill);
     
    % Calculate mean tidal sill speed
    %usill = omega * Aupstream(isill) * a0 / SillData.Xarea(isill);
    %usill = 4 * Aupstream(isill) * a0 / (SillData.Xarea(isill) * Tm2);
    usill = 2 * pi * Aupstream(isill) * a0 / (SillData.Xarea(isill) * Tm2);
    
    % Set flux Richardson number
    Rf = 0.05;
    % if sill flow is super-critical, set Rf = 0.01
    if usill > ci; 
        %disp(['U = ',num2str(usill),' C = ',num2str(ci),' Supercritical flow: Rf = 0.01');
        Rf = 0.01; 
    end
    
    % Calculate work done against buoyancy
    W1 = 0;
    if At > 0
        W1 = Rf * Ej / At;
    end
    
    % Sum work being done
    W_basin = W_basin + W1;

    % Sills generate internal waves propagating in both directions,
    % so calculate seaward energy flux for internal sills.
    if isill < Nsill
        % Set Ht for landward sill
        Ht2 = SillData.Hmax(isill+1);

        % Following lines commented out on 8 November 2011
        % If the surface layer is greater than sill depth, set
        % Ht = Hmean(1)
        %if Hlay(1) > Ht2
        %    Ht2 = Hlay(1);
        %end
        
        % Following lines inserted on 8 November 2011
        if Hlay(1) < Ht2
            Ht2 = Hlay(1);
        end
    
        % Get the mean depth of the water below Ht
        uplayer = floor(Ht2) + 1;
        dnlayer = ceil(Ht2) + 1;
        f = Ht2 - uplayer + 1;
        At_loch = (1-f)*Hypso.A(uplayer) + f*Hypso.A(dnlayer);
        Vt_loch = (1-f)*Hypso.vol(uplayer) + f*Hypso.vol(dnlayer);
        uplayer = floor(Hb) + 1;
        dnlayer = ceil(Hb) + 1;
        f = Hb - uplayer + 1;
        Vb_loch = (1-f)*Hypso.vol(uplayer) + f*Hypso.vol(dnlayer) ...
                              - Vt_loch;
        if Basinarea(isill) > 0
            At = At_loch * Basinarea(isill) / Hypso.A(1);
            Vb = Vb_loch * Basinarea(isill) / Hypso.A(1);    
        else
            % otherwise use whole loch values to get mean depth
            At = At_loch;
            Vb = Vb_loch;
        end
        Hb2 = Vb / At;
    
        % Reduced gravity
        delta_rho = rho(3) - rho(2);
        if Hlay(1) == Ht2
            delta_rho = rho(2) - rho(1);
        end
        %delta_rho = max([delta_rho 0.1]);
        if delta_rho < 0.1
            W(isill) = 0;
        continue
        end
        gprime = delta_rho * Const.g / rho(3);

        % calculate the internal wave phase speed, ci 
        ci = sqrt(gprime * Ht2 * Hb2 / (Ht2 + Hb2));

        % Calculate the energy flux due to internal
        % tide propagation from outer sill
        Ej = 0.5 * rho0 * omega^2 * a0^2 * ci * (Hb2 / (Ht2 + Hb2)) ...
             * Aupstream(isill+1)^2 / SillData.Xarea(isill+1);
    
        % Calculate mean tidal sill speed
        usill = 4 * Aupstream(isill+1) * a0 / ...
            (SillData.Xarea(isill+1) * Tm2);

        % Set flux Richardson number
        Rf = 0.05;
        % if sill flow is super-critical, set Rf = 0.01
        if usill > ci; 
            %disp(['U = ',num2str(usill),' C = ',num2str(ci),' Supercritical flow: Rf = 0.01');
            Rf = 0.01; 
        end
    
        % Calculate the mean rate work done against buoyancy
        W2 = 0;
        if At > 0
            W2 = Rf * Ej / At;
        end
    
        % Sum work being done from all sills in the basin
        W_basin = W_basin + W2;

    end
    
    W(isill) = W_basin;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the rate of density reduction R in the basin. 
    % Calculated based on the total sum of work being done against 
    % buoyancy, assumed to be applied to the water column below the 
    % deeper of the two sills.
    Ht = Ht1;
    Hb = Hb1;
    if isill < Nsill
        Ht = max([Ht1 Ht2]);
        Hb = min([Hb1 Hb2]);
    end

    R(isill) = -Cw * W(isill) / (Const.g * Hb * Hb);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate time scale for renewal events
    % Te = Re / R
    Te(isill) = Re(isill) / abs(R(isill));

    % Convert to days
    Te(isill) = Te(isill) / (3600 * 24);
    
    % Set lower limit
    if Te(isill) < 1; Te(isill) = 1; end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
