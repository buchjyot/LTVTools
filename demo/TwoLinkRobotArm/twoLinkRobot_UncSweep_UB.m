%% Load Data
% This file computes worst-case gain for different uncertainty level and
% compares nominal and robust controller.
load('twoLinkRobot_SpecifyOptions.mat','tvwcopt','tvnopt');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_Hinfsyn.mat');
load('twoLinkRobot_Robsyn.mat','Krob1','Krob2');

% Form Closedloop
CLn0 = lft(evalt(Gunc,Knom.Time),Knom);
CLr1 = lft(evalt(Gunc,Krob1.Time),Krob1);
CLr2 = lft(evalt(Gunc,Krob2.Time),Krob2);

%% Worst-Case Gain Upper Bound Analysis

% Uncertainty Level for Worst-case Gain Analysis
UL = 0:0.1:0.9;
NUL = length(UL);
wcgain0 = zeros(NUL,1);
wcgain1 = zeros(NUL,1);
wcgain2 = zeros(NUL,1);
info0 = cell(NUL,1);
info1 = cell(NUL,1);
info2 = cell(NUL,1);

% Main for loop
parfor i = 1:length(UL)
    tvwcopt1 = tvwcopt;
    tvwcopt1.Display = 'off';
    if i == 1
        [wcg,dWc0,info0{i}] = tvnormb(lft(0,CLn0),NE,tvnopt);wcgain0(i) = wcg(2);
        [wcg,dWc1,info1{i}] = tvnormb(lft(0,CLr1),NE,tvnopt);wcgain1(i) = wcg(2);
        [wcg,dWc2,info2{i}] = tvnormb(lft(0,CLr2),NE,tvnopt);wcgain2(i) = wcg(2);
    else
        tvwcopt1.ULevel = UL(i);
        [wcgain0(i),info0{i}] = tvwcgain(CLn0,Delta,NE,tvwcopt1);
        [wcgain1(i),info1{i}] = tvwcgain(CLr1,Delta,NE,tvwcopt1);
        [wcgain2(i),info2{i}] = tvwcgain(CLr2,Delta,NE,tvwcopt1);
    end
    fprintf(' Completed worst-case gain UB analysis for UL = %.2f\n',UL(i));
end

%% Save Data
save(mfilename,'wcgain0','wcgain1','wcgain2','UL');
return;

%% Plot Data
load('twoLinkRobot_UncSweep_UB.mat');
figure;clf;grid on;box on;hold on;
plot(UL,wcgain0,'-ob',UL,wcgain1,'-rs','LineWidth',2);
legend('$\tilde{T}_0$','$\tilde{T}_{0.8}$',...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
XL = xlim;
YL = ylim;

figure;clf;grid on;box on;hold on;
plot(UL(1),wcgain1(1),'-ob',UL(9),wcgain2(9),'-rs','LineWidth',2);
legend('$\tilde{T}_0$','$\tilde{T}_{0.8}$',...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
xlim(XL);ylim(YL);

figure;clf;grid on;box on;hold on;
plot(UL,wcgain0,'-ob',UL,wcgain2,'-c^',UL,wcgain1,'-rs','LineWidth',2);
legend('$\tilde{T}_0$','$\tilde{T}_{0.4}$','$\tilde{T}_{0.8}$',...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
XL = xlim;
YL = ylim;