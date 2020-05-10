%% Load Data
% This file performs nominal analysis and compares both robust and nominal
% controllers
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Robsyn.mat','Gunc','robinfo','Krob');

%% Nominal Analysis - Nominal Controller vs Robust Controller

% Form Closedloop
CLn = lft(evalt(Gunc,Knom.Time),Knom);
CLr = lft(evalt(Gunc,Krob.Time),Krob);

% Remove Design Weights
CLn = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLn*blkdiag(eye(Nw),inv(Wd),inv(Wn));
CLr = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLr*blkdiag(eye(Nw),inv(Wd),inv(Wn));

% Consider d1,d2 to th1,th2
CLnA = CLn(4:5,2:3);
CLrA = CLr(4:5,2:3);

% Compute worst-case disturbance
[gCLn,dWcCLn] = tvnorm(CLnA,NE,tvnopt);
[gCLr,dWcCLr] = tvnorm(CLrA,NE,tvnopt);

% Normalize worst-case disturbance
dWcCLn = dScl*dWcCLn/tvnorm(dWcCLn);
dWcCLr = dScl*dWcCLr/tvnorm(dWcCLr);

%% Nominal LTV Simulations 
% Response to worst-case disturbance
XCLn = tvlsim(CLnA,dWcCLn,tvspt);
fprintf(' Closed Loop with Nominal Controller :\n');
fprintf(' L2toE Gain Lower Bound = %.3f\n',gCLn(1)*dScl);
fprintf(' LTV Simulation Terminal Euclidean Norm = %.3f\n',norm(tvsubs(XCLn,XCLn.Time(end))));

XCLr = tvlsim(CLrA,-dWcCLr,tvspt);
fprintf(' Closed Loop with Robust Controller :\n');
fprintf(' L2toE Gain Lower Bound = %.3f\n',gCLr(1)*dScl);
fprintf(' LTV Simulation Terminal Euclidean Norm = %.3f\n',norm(tvsubs(XCLr,XCLr.Time(end))));

%% Save Data
save(mfilename,'XCLn','XCLr','gCLn','gCLr','dWcCLn','dWcCLr');