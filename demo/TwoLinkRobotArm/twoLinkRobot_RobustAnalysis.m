%% Load Data
% This file performs robust analysis and compares both robust and nominal
% controllers
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Robsyn.mat','Gunc','Krob1');

%% Robust Analysis - Nominal Controller vs Robust Controller

% Uncertainty Level
tvwcopt.ULevel = DelNorm; % 0.8

% Form Closedloop
CLn = lft(evalt(Gunc,Knom.Time),Knom);
CLr = lft(evalt(Gunc,Krob1.Time),Krob1);

CLn_d_L2 = CLn(1:3,1:3);
CLr_d_L2 = CLr(1:3,1:3);
CLn_d_E  = CLn([1 4:5],1:3);
CLr_d_E  = CLr([1 4:5],1:3);
CLn_d_Total  = CLn(1:5,1:3);
CLr_d_Total  = CLr(1:5,1:3);

CLn_n_L2 = CLn(1:3,[1 4:5]);
CLr_n_L2 = CLr(1:3,[1 4:5]);
CLn_n_E = CLn([1 4:5],[1 4:5]);
CLr_n_E = CLr([1 4:5],[1 4:5]);
CLn_n_Total = CLn(1:5,[1 4:5]);
CLr_n_Total = CLr(1:5,[1 4:5]);

CLn_gd_L2 = CLn(1:3,1:5);
CLr_gd_L2 = CLr(1:3,1:5);
CLn_gd_E = CLn([1 4:5],1:5);
CLr_gd_E = CLr([1 4:5],1:5);
CLn_gd_Total = CLn(1:5,1:5);
CLr_gd_Total = CLr(1:5,1:5);

% Compute worst-case induced gains (Closed loop with nominal controller)
[gCLn_d_L2,dWcCLn_d_L2] = tvwcgain(CLn_d_L2,0,tvwcopt);
[gCLn_d_E,dWcCLn_d_E] = tvwcgain(CLn_d_E,NE,tvwcopt);
[gCLn_d_Total,dWcCLn_d_Total] = tvwcgain(CLn_d_Total,NE,tvwcopt);

[gCLn_n_L2,dWcCLn_n_L2] = tvwcgain(CLn_n_L2,0,tvwcopt);
[gCLn_n_E,dWcCLn_n_E] = tvwcgain(CLn_n_E,NE,tvwcopt);
[gCLn_n_Total,dWcCLn_n_Total] = tvwcgain(CLn_n_Total,NE,tvwcopt);

[gCLn_gd_L2,dWcCLn_gd_L2] = tvwcgain(CLn_gd_L2,0,tvwcopt);
[gCLn_gd_E,dWcCLn_gd_E] = tvwcgain(CLn_gd_E,NE,tvwcopt);
[gCLn_gd_Total,dWcCLn_gd_Total] = tvwcgain(CLn_gd_Total,NE,tvwcopt);

% Compute worst-case induced gains (Closed loop with robust controller)
[gCLr_d_L2,dWcCLr_d_L2] = tvwcgain(CLr_d_L2,0,tvwcopt);
[gCLr_d_E,dWcCLr_d_E] = tvwcgain(CLr_d_E,NE,tvwcopt);
[gCLr_d_Total,dWcCLr_d_Total] = tvwcgain(CLr_d_Total,NE,tvwcopt);

[gCLr_n_L2,dWcCLr_n_L2] = tvwcgain(CLr_n_L2,0,tvwcopt);
[gCLr_n_E,dWcCLr_n_E] = tvwcgain(CLr_n_E,NE,tvwcopt);
[gCLr_n_Total,dWcCLr_n_Total] = tvwcgain(CLr_n_Total,NE,tvwcopt);

[gCLr_gd_L2,dWcCLr_gd_L2] = tvwcgain(CLr_gd_L2,0,tvwcopt);
[gCLr_gd_E,dWcCLr_gd_E] = tvwcgain(CLr_gd_E,NE,tvwcopt);
[gCLr_gd_Total,dWcCLr_gd_Total] = tvwcgain(CLr_gd_Total,NE,tvwcopt);

%% Save
save(mfilename);