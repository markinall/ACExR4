%% Script to set up your environment and run the ACExR model
% NOTES: MEI March 2021
% All m-files are in the "code" directory
% All loch topog data are csv files in the "Catalogue" directory
% All forcing data are csv files in the "in" directory
% 
% ACExR.m is the main model run script/function, it calls various functions:
% 
% 1) ACconfigure.m - this sets the switches to configure the model. RUN
% LENGTH IN DAYS IS SET HERE (L40).
% 
% 2) ReadCatalogue.m - this reads Read topographic data from sea loch catalogue
% 
% 3) Hypsography.m - Derive a hypsographic function for the loch
% 
% 4) ReadDBForcing.m - reads forcing data
% 
% 5) Initialise.m - Initialises model parameters from topographic and boundary forcing data.
% 
% 6) CalcE.m - runs the physics functions of the model
%
% 7) Output as structured arrays in a mat file called results.mat in "code"
% directory (not the best choice of location!)

% House keeping
close all
clear
clc

% set the paths (only need to do this once per matlab session, but included
% here for completness
script_to_set_the_path

start_folder = pwd; % save start location
% cd to code sub-dir
if ismac
    cd ./code
else
    cd .\code
end
% now run model for loch Creran
ACExR('Creran')
clear cd
cd(start_folder) % return to start location




