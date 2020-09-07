%% Load Data
% This file performs nominal analysis and compares both robust and nominal
% controllers
load('TwoStateEx.mat');
load('TwoStateEx_NominalSynthesis.mat');
load('TwoStateEx_RobustSynthesis.mat');

%% Nominal Analysis - Nominal Controller vs Robust Controller
% Remove Design Weights
CLn = blkdiag(eye(Nv),inv(Weu),inv(Wex))*CLn*blkdiag(eye(Nw),inv(Wd),inv(Wn));
CLr = blkdiag(eye(Nv),inv(Weu),inv(Wex))*CLr*blkdiag(eye(Nw),inv(Wd),inv(Wn));

% Consider d1 to eE
CLnA = CLn(3:4,2);
CLrA = CLr(3:4,2);

% Compute worst-case disturbance
[gCLn,dWcCLn] = tvnorm(CLnA,NE,tvnopt);
[gCLr,dWcCLr] = tvnorm(CLrA,NE,tvnopt);

% Normalize worst-case disturbance
dScl = 1;
dWcCLn = dScl*dWcCLn/tvnorm(dWcCLn);
dWcCLr = dScl*dWcCLr/tvnorm(dWcCLr);

% Response to worst-case disturbance
XCLn = tvlsim(CLnA,dWcCLn,tvopt);
fprintf(' Closed Loop with Nominal Controller :\n');
fprintf(' L2toE Gain Lower Bound = %.3f\n',gCLn(1)*dScl);
fprintf(' LTV Simulation Terminal Euclidean Norm = %.3f\n',norm(tvsubs(XCLn,Tf)));

XCLr = tvlsim(CLrA,dWcCLr,tvopt);
fprintf(' Closed Loop with Robust Controller :\n');
fprintf(' L2toE Gain Lower Bound = %.3f\n',gCLr(1)*dScl);
fprintf(' LTV Simulation Terminal Euclidean Norm = %.3f\n',norm(tvsubs(XCLr,Tf)));

%% Save Data
save(mfilename,'gCLn','gCLn','gCLr','dScl','XCLn','XCLr','T0','Tf','beta');