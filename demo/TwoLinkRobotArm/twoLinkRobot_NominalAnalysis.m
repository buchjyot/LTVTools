%% Load Data
% This file performs nominal analysis and compares both robust and nominal
% controllers
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Robsyn.mat','Gunc','Krob1');

%% Nominal Analysis - Nominal Controller vs Robust Controller

% Form Closedloop
CLn = lft(evalt(Gunc,Knom.Time),Knom);
CLr = lft(evalt(Gunc,Krob1.Time),Krob1);

CLn_d_L2 = CLn(2:3,2:3);
CLr_d_L2 = CLr(2:3,2:3);
CLn_d_E = CLn(4:5,2:3);
CLr_d_E = CLr(4:5,2:3);
CLn_d_Total = CLn(2:5,2:3);
CLr_d_Total = CLr(2:5,2:3);

CLn_n_L2 = CLn(2:3,4:5);
CLr_n_L2 = CLr(2:3,4:5);
CLn_n_E = CLn(4:5,4:5);
CLr_n_E = CLr(4:5,4:5);
CLn_n_Total = CLn(2:5,4:5);
CLr_n_Total = CLr(2:5,4:5);

CLn_gd_L2 = CLn(2:3,2:5);
CLr_gd_L2 = CLr(2:3,2:5);
CLn_gd_E = CLn(4:5,2:5);
CLr_gd_E = CLr(4:5,2:5);
CLn_gd_Total = CLn(2:5,2:5);
CLr_gd_Total = CLr(2:5,2:5);

% Compute induced gains (Closed loop with nominal controller)
[gCLn_d_L2,dWcCLn_d_L2] = tvnormb(CLn_d_L2,0,tvnopt);
[gCLn_d_E,dWcCLn_d_E] = tvnormb(CLn_d_E,NE,tvnopt);
[gCLn_d_Total,dWcCLn_d_Total] = tvnormb(CLn_d_Total,NE,tvnopt);

[gCLn_n_L2,dWcCLn_n_L2] = tvnormb(CLn_n_L2,0,tvnopt);
[gCLn_n_E,dWcCLn_n_E] = tvnormb(CLn_n_E,NE,tvnopt);
[gCLn_n_Total,dWcCLn_n_Total] = tvnormb(CLn_n_Total,NE,tvnopt);

[gCLn_gd_L2,dWcCLn_gd_L2] = tvnormb(CLn_gd_L2,0,tvnopt);
[gCLn_gd_E,dWcCLn_gd_E] = tvnormb(CLn_gd_E,NE,tvnopt);
[gCLn_gd_Total,dWcCLn_gd_Total] = tvnormb(CLn_gd_Total,NE,tvnopt);

% Compute induced gains (Closed loop with robust controller)
[gCLr_d_L2,dWcCLr_d_L2] = tvnormb(CLr_d_L2,0,tvnopt);
[gCLr_d_E,dWcCLr_d_E] = tvnormb(CLr_d_E,NE,tvnopt);
[gCLr_d_Total,dWcCLr_d_Total] = tvnormb(CLr_d_Total,NE,tvnopt);

[gCLr_n_L2,dWcCLr_n_L2] = tvnormb(CLr_n_L2,0,tvnopt);
[gCLr_n_E,dWcCLr_n_E] = tvnormb(CLr_n_E,NE,tvnopt);
[gCLr_n_Total,dWcCLr_n_Total] = tvnormb(CLr_n_Total,NE,tvnopt);

[gCLr_gd_L2,dWcCLr_gd_L2] = tvnormb(CLr_gd_L2,0,tvnopt);
[gCLr_gd_E,dWcCLr_gd_E] = tvnormb(CLr_gd_E,NE,tvnopt);
[gCLr_gd_Total,dWcCLr_gd_Total] = tvnormb(CLr_gd_Total,NE,tvnopt);

%% Save
save(mfilename);