%% Robust Analysis of Doyle's Example

% Load Workspace
clear;load('DoyleLQGExampleFH.mat');

%% LQR Robust Analysis
% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
fprintf('==================================================\n');
fprintf('Robust Closed-Loop Analysis [LQR] \n')
fprintf('==================================================\n');

% w to v gain
tvnLQR = tvnorm(Tlqrunc(1,1),tvnopt) %#ok<*NOPTS>

%% LQG Robust Analysis
% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
fprintf('==================================================\n');
fprintf('Robust Closed-Loop Analysis [LQG] \n')
fprintf('==================================================\n');

% w to v gain
tvnLQG1 = tvnorm(Tlqgunc(1,1),tvnopt)

% Noise to output gain
tvnLQG2 = tvnorm(Tlqgunc(2,4),tvnopt)
