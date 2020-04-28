%% twoLinkRobot_HinfDesign.m
% This file performs sysic to do robust synthesis

%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat');

%% Signal Dimentions
% Number of outputs penalized in Euclidean sense
NE = 2;

% Number of measurement available for feedback
Ny = 2;

%% Signal Weights
% Disturbance weight
WdScl = 0.1;
Wd = ss(WdScl*eye(Nu));

% Weights on Euclidean Outputs
WEScl = 1;
WE = ss(WEScl*eye(NE));

% Weights on control effort
WuScl = 0.5;
Wu = ss(WuScl*eye(Nu));

% Noise weight
WnScl = 0.01;
Wn = ss(WnScl*eye(Ny));

% Uncertainty Channel Input Scalling (assumed to be Identity)
WvScl = 1;
Wv = ss(WvScl*eye(Nv));

% Uncertainty Channel Output Scalling (assumed to be Identity)
WwScl = 1;
Ww = ss(WwScl*eye(Nw));

%% State-Feedback Problem
% Assumes that we have an access to all the states
systemnames = 'G Wd Wu WE'; %#ok<*NASGU>
inputvar = '[d(2); u(2)]';
outputvar = '[Wu; WE]';
input_to_G = '[u+Wd]';
input_to_Wd = '[d]';
input_to_Wu = '[u]';
input_to_WE = '[G(1:2)]';
cleanupsysic = 'yes';
Gsfb = sysic;

%% Uncertain Measurement-Feedback Problem
% Input Uncertainity enters in the second control channel
systemnames = 'G Wv Ww Wd Wn Wu WE'; %#ok<*NASGU>
inputvar = '[w; d(2); n(2); u(2)]';
outputvar = '[Wv; Wu; WE; G(1:2)+Wn]';
input_to_G = '[u(1)+Wd(1);u(2)+Wd(2)+Ww]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
input_to_Wu = '[u]';
input_to_WE = '[G(1:2)]';
input_to_Wv = '[Wu(2)+Wd(2)]';
input_to_Ww = '[w]';
cleanupsysic = 'yes';
Gunc = sysic;

% Plant to be used for nominal synthesis
Gnom = lft(0,Gunc);

%% Save Data
save(mfilename,'Gsfb','Gnom','Gunc','xhat0','Wd','Wn','Wu','Wv','Ww','WE','NE','Ny');