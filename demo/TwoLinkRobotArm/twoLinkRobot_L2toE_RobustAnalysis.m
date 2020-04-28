%% Load Data
% This file performs robust analysis and compares both robust and nominal
% controllers
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Robsyn.mat','Gunc','robinfo','Krob');

%% Robust Analysis - Nominal Controller vs Robust Controller

% Form Closedloop
CLn = lft(evalt(Gunc,Knom.Time),Knom);
CLr = lft(evalt(Gunc,Krob.Time),Krob);

% Remove Design Weights
CLn = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLn*blkdiag(eye(Nw),inv(Wd),inv(Wn));
CLr = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLr*blkdiag(eye(Nw),inv(Wd),inv(Wn));

% Consider d1,d2 to th1,th2
CLnA = CLn([1 4:5],1:3);
CLrA = CLr([1 4:5],1:3);

% Uncertainty Norm Bound
tvwcopt.ULevel = DelNorm;

% Worst-case Gain
[wcgUB1,info1] = tvwcgain(CLnA,Delta,NE,tvwcopt);
[wcgUB2,info2] = tvwcgain(CLrA,Delta,NE,tvwcopt);

%% Sample Uncertainties
NSamples = 100;
rng('shuffle');
UncSample = ltiusample(DelNorm,NSamples);
rng('default');
UncSample = [UncSample;DelNorm;-DelNorm];

%% Find Worst-Case Unceratinty
N = length(UncSample);
g1 = zeros(N,1);
g2 = zeros(N,1);
dWc1 = cell(N,1);
dWc2 = cell(N,1);
parfor i = 1:N
    [gE1,dWc1{i}] = tvnorm(lft(UncSample{i},CLnA),NE,tvnopt);
    g1(i) = gE1(1);
    
    [gE2,dWc2{i}] = tvnorm(lft(UncSample{i},CLrA),NE,tvnopt);
    g2(i) = gE2(1);
    
    fprintf(' wcg1:%.3f, wcg2:%.3f\n',g1(i),g2(i));
end

% Find worst-case uncertainty
[wcgLB1,id1] = max(g1);
[wcgLB2,id2] = max(g2);
Deltawc1 = UncSample{id1};
Deltawc2 = UncSample{id2};

% Display results
fprintf(' Nominal Controller Closed Loop Worst-Case Gain: [%.4f,%.4f]\n',wcgLB1,wcgUB1);
fprintf(' Worst-Case Delta:\n');
Deltawc1 %#ok<*NOPTS>

fprintf(' Robust Controller Closed Loop Worst-Case Gain: [%.4f,%.4f]\n',wcgLB2,wcgUB2);
fprintf(' Worst-Case Delta:\n');
Deltawc2

%% LTV Simulations with Worst-Case Uncertainty and Worst-Case Disturbance
% Normalize worst-case disturbance
dWc1s = dScl*dWc1{id1}/tvnorm(dWc1{id1});
dWc2s = dScl*dWc2{id2}/tvnorm(dWc2{id2});

% Response to worst-case disturbance & worst-case uncertainty
XCLn = tvlsim(lft(Deltawc1,CLnA),dWc1s,tvopt);
XCLr = tvlsim(lft(Deltawc2,CLrA),-dWc2s,tvopt);

%% Save Data
save(mfilename,'Deltawc1','Deltawc2','wcgUB1','wcgUB2','wcgLB1','wcgLB2','XCLn','XCLr','dWc1s','dWc2s');