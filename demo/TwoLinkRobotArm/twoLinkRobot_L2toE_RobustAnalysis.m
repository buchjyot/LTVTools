%% Load Data
% This file performs robust analysis and compares both robust and nominal
% controllers
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Robsyn.mat','Gunc','Krob1');

%% Robust Analysis - Nominal Controller vs Robust Controller

% Form Closedloop
CLn0 = lft(evalt(Gunc,Knom.Time),Knom);
CLr1 = lft(evalt(Gunc,Krob1.Time),Krob1);

% Remove Design Weights
% CLn0 = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLn0*blkdiag(eye(Nw),inv(Wd),inv(Wn));
% CLr1 = blkdiag(eye(Nv),inv(Wu),inv(WE))*CLr1*blkdiag(eye(Nw),inv(Wd),inv(Wn));

% Consider d1,d2 to th1,th2
CLn0A = CLn0([1 4:5],:);
CLr1A = CLr1([1 4:5],:);

% Worst-case Gain
tvwcopt.ULevel = DelNorm; % 0.8
[wcgUB0,info0] = tvwcgain(CLn0A,Delta,NE,tvwcopt);
[wcgUB1,info1] = tvwcgain(CLr1A,Delta,NE,tvwcopt);

%% Sample Uncertainties
NSamples = 100;
UncSample = ltiusample(DelNorm,NSamples);
UncSample = [UncSample;DelNorm;-DelNorm];

% Find Worst-Case Unceratinty
N = length(UncSample);
g0 = zeros(N,1);
g1 = zeros(N,1);
dWc0 = cell(N,1);
dWc1 = cell(N,1);
parfor i = 1:N
    [gE0,dWc0{i}] = tvnormb(tvss(lft(UncSample{i},CLn0A.Data),CLn0A.Time),NE,tvnopt);
    g0(i) = gE0(1);
    
    [gE1,dWc1{i}] = tvnormb(tvss(lft(UncSample{i},CLr1A.Data),CLr1A.Time),NE,tvnopt);
    g1(i) = gE1(1);
    
    fprintf(' wcg0:%.3f, wcg1:%.3f\n',g0(i),g1(i));
end

% Find worst-case uncertainty
[wcgLB0,id0] = max(g0);
[wcgLB1,id1] = max(g1);
Deltawc0 = UncSample{id0};
Deltawc1 = UncSample{id1};

% Display results
fprintf(' Nominal Controller Closed Loop Worst-Case Gain: [%.4f,%.4f]\n',wcgLB0,wcgUB0);
fprintf(' Worst-Case Delta:\n');
Deltawc0 %#ok<*NOPTS>

fprintf(' Robust Controller Closed Loop Worst-Case Gain: [%.4f,%.4f]\n',wcgLB1,wcgUB1);
fprintf(' Worst-Case Delta:\n');
Deltawc1

%% LTV Simulations with Worst-Case Uncertainty and Worst-Case Disturbance
% Normalize worst-case disturbance
dWc0s = dScl*dWc0{id0}/tvnorm(dWc0{id0});
dWc1s = dScl*dWc1{id1}/tvnorm(dWc1{id1});

% Response to worst-case disturbance & worst-case uncertainty
XCLn = tvlsim(tvss(lft(Deltawc0,CLn0A.Data),CLn0A.Time),dWc0s,tvopt);
XCLr = tvlsim(tvss(lft(Deltawc1,CLr1A.Data),CLr1A.Time),-dWc1s,tvopt);

%% Save Data
save(mfilename,'Deltawc0','Deltawc1','wcgUB0','wcgUB1','wcgLB0','wcgLB1','XCLn','XCLr','dWc0s','dWc1s');