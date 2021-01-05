%% DEMO2
% This example computes the worst-case norm of a matrix using mixed mu
% power iteration code in MATLAB.

% Create random matrix
rng(1);
M = rand(4);

% Uncertainty
Del = ucomplexm('Del',zeros(2));

% Interconnection
P = lft(Del,M);

% Worst-case norm
[wcn,wcu] = wcnorm(P);

% NOTE: the power iteration code in shipping version of MATLAB is provided
% in edit(fullfile(matlabroot,'toolbox','robust','rctutil','private','mmupiter.m'))