%% Load Data
% This file computes worst-case gain for different uncertainty level and
% compares nominal and robust controller.
if false
    load('twoLinkRobot_SpecifyOptions.mat','tvwcopt','tvnopt');
    load('twoLinkRobot_HinfDesign.mat');
    load('twoLinkRobot_Hinfsyn.mat');
    load('twoLinkRobot_Robsyn.mat','Krob1','Krob2');
    
    % Form Closedloop
    CLn0 = lft(evalt(Gunc,Knom.Time),Knom);
    CLr1 = lft(evalt(Gunc,Krob1.Time),Krob1);
    CLr2 = lft(evalt(Gunc,Krob2.Time),Krob2);
    save(mfilename);
else
    load('twoLinkRobot_UncSweep_LB_Data.mat');
end

%% Nominal Analysis
% Lower bound at nominal
if false
    wcg = tvnormb(tvss(lft(0,CLn0.Data),CLn0.Time),NE,tvnopt);wcg00 = wcg(1);
    wcg = tvnormb(tvss(lft(0,CLr1.Data),CLr1.Time),NE,tvnopt);wcg11 = wcg(1);
    wcg = tvnormb(tvss(lft(0,CLr2.Data),CLr2.Time),NE,tvnopt);wcg22 = wcg(1);
    save([mfilename '_0']);
end

%% Robust Analysis

if true
    % Uncertainty Level for Worst-case Gain Analysis
    UL = 0.9;
    NUL = length(UL);
    
    % Sample Uncertainties
    NSamples = 100;
    UncSample = ltiusample(UL,NSamples);
    
    %% Main for loop
    wcg0 = zeros(NSamples,1);
    wcg1 = zeros(NSamples,1);
    wcg2 = zeros(NSamples,1);
    parpool(4);
    parfor i = 1:NSamples
        fprintf(' j=%d\n',i);
        Delta = UncSample{i};
        
        temp0 = tvss(lft(Delta,CLn0.Data),CLn0.Time);
        wcg = tvnormb(temp0,NE,tvnopt);wcg0(i) = wcg(1);
        
        temp1 = tvss(lft(Delta,CLr1.Data),CLr1.Time);
        wcg = tvnormb(temp1,NE,tvnopt);wcg1(i) = wcg(1);
        
        temp2 = tvss(lft(Delta,CLr2.Data),CLr2.Time);
        wcg = tvnormb(temp2,NE,tvnopt);wcg2(i) = wcg(1);
    end
    
    [wcgLB0,id0] = max(wcg0); Deltawc0 = UncSample{id0};
    [wcgLB1,id1] = max(wcg1); Deltawc1 = UncSample{id1};
    [wcgLB2,id2] = max(wcg2); Deltawc2 = UncSample{id2};
    fprintf(' Completed worst-case gain LB analysis for UL = %.2f\n',UL);
    
    %% Save Data
    save([mfilename '_' num2str(UL*10)],'UL','wcg0','wcg1','wcg2','wcgLB0','wcgLB1','wcgLB2','Deltawc0','Deltawc1','Deltawc2');
end
return;

%% Plot Lower Bounds
UncLevel = 0:0.1:0.9;
NUL = length(UncLevel);

% Memory allocation
wcgLB0All = zeros(NUL,1);
wcgLB1All = zeros(NUL,1);
wcgLB2All = zeros(NUL,1);
wcDelta0  = cell(NUL,1);
wcDelta1  = cell(NUL,1);
wcDelta2  = cell(NUL,1);

% For loop
for i = 1:NUL
    load(['twoLinkRobot_UncSweep_LB_' num2str(UncLevel(i)*10) '.mat'])
    if i==1
        wcgLB0All(i) = wcg00;wcDelta0{i} = 0;
        wcgLB1All(i) = wcg11;wcDelta1{i} = 0;
        wcgLB2All(i) = wcg22;wcDelta2{i} = 0;
    else
        wcgLB0All(i) = wcgLB0;wcDelta0{i} = Deltawc0;
        wcgLB1All(i) = wcgLB1;wcDelta1{i} = Deltawc1;
        wcgLB2All(i) = wcgLB2;wcDelta2{i} = Deltawc2;
    end
end

% Save Data
UL = UncLevel;
save('twoLinkRobot_UncSweep_LB.mat','UL','wcgLB0All','wcgLB1All','wcgLB2All',...
    'wcDelta0','wcDelta1','wcDelta2');

%% Plot Lower Bounds
load('twoLinkRobot_UncSweep_LB.mat');

figure;clf;grid on;box on;hold on;
plot(UL,wcgLB0All,'--ob',UL,wcgLB1All,'--c^',UL,wcgLB2All,'--rs','LineWidth',2);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
XL = xlim;
YL = ylim;

%% Plot upper bound and lower bounds together
load('twoLinkRobot_UncSweep_LB.mat');
load('twoLinkRobot_UncSweep_UB.mat');

figure;clf;grid on;box on;hold on;
plot(UL,wcgain0,'-ob',UL,wcgain2,'-c^',UL,wcgain1,'-rs','LineWidth',2);
plot(UL,wcgLB0All,'--ob',UL,wcgLB2All,'--c^',UL,wcgLB1All,'--rs','LineWidth',2);
legend('$\tilde{T}_0$','$\tilde{T}_{0.4}$','$\tilde{T}_{0.8}$',...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
XL = xlim;
YL = ylim;