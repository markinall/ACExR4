function ACExR(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ACExR(varargin)
%
% Parameterisation of physical Exchange Rates in
% Scottish fjords for Assimilative Capacity modelling.
%
% Usage:    E contains exchange rates (units s-1) of each layer
%               and inter-layer exchange. 
%           P contains parameter values e.g. layer volumes, 
%               thicknesses, salinity, temperature etc.
%           Specify the loch to model, and the corresponding databases by
%               editting ACconfigure.m
%
% Phil Gillibrand
% 10/3/2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variables
global LochData SillData Hypso Bdata Const D E Param

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise arrays
LochData = [];
SillData = [];
Hypso = [];
Bdata = [];
Const = [];
D = [];
E = [];
Param = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read model configuration data. If a name is supplied use it, otherwise
% look to the configuration file for a name.
nargin = length(varargin);
if nargin > 0
    LochData.Name = varargin{1};
    ACconfigure(LochData.Name);
else
    ACconfigure;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read topographic data from sea loch catalogue
% If riverine nutrient loadings are supplied as arguments (tonnes/year),
% overwrite the values loaded from the database (converting to kg/s).
ReadCatalogue;
if nargin > 1
    LochData.QC1 = str2double(varargin{2}) * 1000 /(365 * 86400);
end
if nargin > 2
    LochData.QC2 = str2double(varargin{3}) * 1000 /(365 * 86400);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derive a hypsographic function for the loch
Hypsography;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read forcing data
ReadDBForcing2;
% save('ABC',Bdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise model parameters
Initialise;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run model for annual cycle
try
    CalcE;
catch
    %if Param.Error > 0
        %disp(['ERROR: Run stopped. Error code = ',num2str(Param.Error)])
    % end
    clear D;
    dt = Const.deltaT;
    Initialise;
    Const.deltaT = dt / 2;
    CalcE;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
end

