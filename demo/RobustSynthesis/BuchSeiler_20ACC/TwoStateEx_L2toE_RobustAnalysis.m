%% Load Data
% This file performs robust analysis and compares both robust and nominal
% controllers
load('TwoStateEx.mat');
load('TwoStateEx_NominalSynthesis.mat');
load('TwoStateEx_RobustSynthesis.mat');
load('TwoStateEx_SampleLTIUnc');

%% Robust Analysis - Nominal Controller vs Robust Controller
% Remove Design Weights
CLn = blkdiag(eye(Nv),inv(Weu),inv(Wex))*CLn*blkdiag(eye(Nw),inv(Wd),inv(Wn));
CLr = blkdiag(eye(Nv),inv(Weu),inv(Wex))*CLr*blkdiag(eye(Nw),inv(Wd),inv(Wn));

% Get rid of the design weights but retain uncertainty channels
CLnA = CLn([1 3 4],[1 2]);
CLrA = CLr([1 3 4],[1 2]);

% Worst-case gain
tvwcopt.ULevel = beta;
[wcgUB1,info1] = tvwcgain(CLnA,Delta,NE,tvwcopt);
[wcgUB2,info2] = tvwcgain(CLrA,Delta,NE,tvwcopt);

%% Sample across uncertainties
NSample = length(UncSample);
g1 = zeros(NSample,1);
g2 = zeros(NSample,1);
dWc1 = cell(NSample,1);
dWc2 = cell(NSample,1);
parfor i = 1:NSample
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

% Normalize worst-case disturbance
dScl = 1;
dWc1s = dScl*dWc1{id1}/tvnorm(dWc1{id1});
dWc2s = dScl*dWc2{id2}/tvnorm(dWc2{id2});

% Response to worst-case disturbance & worst-case uncertainty
XCLn = tvlsim(lft(Deltawc1,CLnA),dWc1s,tvopt);
XCLr = tvlsim(lft(Deltawc2,CLrA),-dWc2s,tvopt);

%% Save Data
save(mfilename,'Deltawc1','Deltawc2','T0','Tf','wcgUB1','wcgUB2','wcgLB1','wcgLB2','XCLn','XCLr','dWc1s','dWc2s','dScl');