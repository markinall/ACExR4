% a_script_to_set_the_path
% 
% run this script once at the commencement of a session
% PT, 10 April 2011, 11 May 2011
% MI, March 2021 for mac os

% display start information in command window
fprintf('\n%s\n', '-----------------------------------------');
fprintf('%s\n', 'Script to set path for ACExR started at: ');
fprintf('%s\n', ['            ' datestr(now)]);

%% following line must be added for computer/OS in use
 cd = './';
% Add sub-directories for the duration of the session
% this is based on code from cdsarfac made by EP in 2006
% and should work on any platform
actpath=cd; % get current path
rec=genpath(actpath); % returns a path string formed by 
% recursively adding all the directories below the current directory.
addpath(rec); % use this string to add all to current search path

fprintf('%s\n', 'Current directory is: ');
fprintf('%s\n', pwd);
fprintf('%s\n', '-- all subdirectories have been added to path');
fprintf('%s\n\n', '------------------------------------------');