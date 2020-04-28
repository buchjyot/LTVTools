%% Load Data
% This file computes worst-case gain for different uncertainty level and
% compares nominal and robust controller.
load('twoLinkRobot_SpecifyOptions.mat','tvwcopt','tvnopt');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_Robsyn.mat');

% Form Closedloop
CLn = lft(evalt(Gunc,Knom.Time),Knom);
CLr = lft(evalt(Gunc,Krob.Time),Krob);

%% Worst-Case Gain Analysis

% Uncertainty Level for Worst-case Gain Analysis
UL = 0:0.1:0.9;
NUL = length(UL);
wcgain1 = zeros(NUL,1);
wcgain2 = zeros(NUL,1);
info1 = cell(NUL,1);
info2 = cell(NUL,1);

% Main for loop
parfor i = 1:length(UL)
    tvwcopt1 = tvwcopt;
    tvwcopt1.Display = 'off';
    if i == 1
        [wcg,dWc1,info1{i}] = tvnorm(lft(0,CLn),NE,tvnopt);wcgain1(i) = wcg(2);
        [wcg,dWc2,info2{i}] = tvnorm(lft(0,CLr),NE,tvnopt);wcgain2(i) = wcg(2);
    else
        tvwcopt1.ULevel = UL(i);
        [wcgain1(i),info1{i}] = tvwcgain(CLn,Delta,NE,tvwcopt1);
        [wcgain2(i),info2{i}] = tvwcgain(CLr,Delta,NE,tvwcopt1);
    end
    fprintf(' Completed worst-case gain analysis for UL = %.2f\n',UL(i));
end

%% Save Data
save(mfilename,'wcgain1','wcgain2','UL');