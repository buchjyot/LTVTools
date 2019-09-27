%% Add LTVTools to MATLAB path

% Add Path to source directories
addpath(ltvroot);
addpath(genpath(fullfile(ltvroot,'toolbox')));
addpath(genpath(fullfile(ltvroot,'test')));
addpath(genpath(fullfile(ltvroot,'demo')));
addpath(genpath(fullfile(ltvroot,'doc')));

% sl_refresh_customizations for ltvlib simulink library
sl_refresh_customizations

% rehash toolbox
rehash toolboxcache;
rehash toolboxreset;