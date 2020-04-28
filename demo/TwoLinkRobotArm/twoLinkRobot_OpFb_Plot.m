%% PLOT LQR and OpFb Robust Gains
% This file plots the robust L2-to-Euclidean gains for LQR and Output
% feedback interconnections. First, control design with LQR, OpFb(Tau =
% 0.1), OpFb(tau = 0.05), OpFb(tau = 0.01) was performed. Then worst-case
% gain was computed as uncertainty level (beta) is varied.
% Where, ||\delta|| <= beta
%
% NOTE: "twoLinkRobot_OpFb" file does one such analysis at beta = 0.1.

% Load the MAT file
load('RobustL2toEGains_LQR_OpFb.mat')

% Define Colors
color1 = 'r';
color2 = 'g';
color3 = 'b';
color4 = 'c';

% Disturbances are scaled to have norm ||d|| = dL2norm
dL2norm = 5;

% Transperency
alphaSetting = 1;

%% Compare Robust Gains
figure;
yLim = 0:0.1:3;

p1 = plot(DelNormArray_OpFb_tau1,gCL_OpFb_tau1,'-*','Color',color1,'LineWidth',2);hold on;
plot(0.2*ones(1,numel(yLim)),yLim,'--','Color',color1,'LineWidth',2);

p2 = plot(DelNormArray_OpFb_tau2,gCL_OpFb_tau2,'-^','Color',color2,'LineWidth',2);hold on;
plot(0.31*ones(1,numel(yLim)),yLim,'--','Color',color2,'LineWidth',2);

p3 = plot(DelNormArray_OpFb_tau3,gCL_OpFb_tau3,'-+','Color',color3,'LineWidth',2);hold on;
plot(0.504*ones(1,numel(yLim)),yLim,'--','Color',color3,'LineWidth',2);

p4 = plot(DelNormArray_LQR,gCL_LQR,'-o','Color',color4,'LineWidth',2);
plot(0.9411*ones(1,numel(yLim)),yLim,'--','Color',color4,'LineWidth',2);

set(gca,'YLim',[0 0.4]);

% Annotations
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst Case Gain (\gamma_{wc})','FontSize',14)
title(sprintf('Worst-Case Robust Gain Plot'))
legend([p1 p2 p3 p4],'\tau = 0.1','\tau = 0.05',...
    '\tau = 0.02','LQR','FontSize',14);
box on;grid on;
alpha(alphaSetting);

%% Figure 2
figure;
yLim = 0:0.1:3;

p1 = plot(DelNormArray_OpFb_tau1,gCL_OpFb_tau1,'-*','Color',color1,'LineWidth',2);hold on;

p2 = plot(DelNormArray_OpFb_tau2,gCL_OpFb_tau2,'-^','Color',color2,'LineWidth',2);hold on;

p3 = plot(DelNormArray_OpFb_tau3,gCL_OpFb_tau3,'-+','Color',color3,'LineWidth',2);hold on;

p4 = plot(DelNormArray_LQR,gCL_LQR,'-o','Color',color4,'LineWidth',2);

set(gca,'YLim',[0 0.4]);

% Annotations
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst Case Gain (\gamma_{wc})','FontSize',14)
title(sprintf('Worst-Case Robust Gain Plot'))
legend([p1 p2 p3 p4],'\tau = 0.1','\tau = 0.05',...
    '\tau = 0.02','LQR','FontSize',14);
set(gca,'XLim',[0 0.2]);
box on;grid on;
alpha(alphaSetting);

%% Plot 1 - LQR vs Nominal Trajectory
% Refresh Workspace
load('twoLinkRobot_BuildLTVModel.mat','etabar','etaf','eta0');
cBall = tvsubs( etabar(1:2), etabar.Time(end));
NBall = 50;
aBall = linspace(0,2*pi,NBall);

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
gCL = gCL_LQR(DelNormArray_LQR == 0.1);
rBall = gCL*dL2norm;
th1Ball1 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball1 = cBall(2) + rBall*sin(aBall);

figure;clf;
patch(th1Ball1,th2Ball1,color4);
hold on;
plot(etabar(1)+2*pi,etabar(2),'k',eta0(1)+2*pi, eta0(2),'ko',...
    etaf(1)+2*pi, etaf(2),'kx','LineWidth',2,'MarkerSize',6);
axis equal;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;box on;grid on;
set(gca,'xlim',[-0.0612    5.0826]);
set(gca,'ylim',[-3.5532    0.5038]);
alpha(alphaSetting);

%% Plot 2 - LQR vs tau = 0.02 disk

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
gCL = gCL_OpFb_tau3(DelNormArray_OpFb_tau3 == 0.1);
rBall = gCL*dL2norm;
th1Ball2 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball2 = cBall(2) + rBall*sin(aBall);

figure;clf;
patch(th1Ball2,th2Ball2,color3);
hold on;patch(th1Ball1,th2Ball1,color4);
plot(etabar(1)+2*pi,etabar(2),'k',eta0(1)+2*pi, eta0(2),'ko',...
    etaf(1)+2*pi, etaf(2),'kx','LineWidth',2,'MarkerSize',6);
axis equal;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;box on;grid on;
set(gca,'xlim',[-0.0612    5.0826]);
set(gca,'ylim',[-3.5532    0.5038]);
alpha(alphaSetting);

%% Plot 3 - LQR vs tau = 0.05 disk

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
gCL = gCL_OpFb_tau2(DelNormArray_OpFb_tau2 == 0.1);
rBall = gCL*dL2norm;
th1Ball3 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball3 = cBall(2) + rBall*sin(aBall);

figure;clf;
patch(th1Ball3,th2Ball3,color2);
hold on;patch(th1Ball1,th2Ball1,color4);
plot(etabar(1)+2*pi,etabar(2),'k',eta0(1)+2*pi, eta0(2),'ko',...
    etaf(1)+2*pi, etaf(2),'kx','LineWidth',2,'MarkerSize',6);
axis equal;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;box on;grid on
set(gca,'xlim',[-0.0612    5.0826]);
set(gca,'ylim',[-3.5532    0.5038]);
alpha(alphaSetting);

%% Plot 4 - LQR vs tau = 0.02 disk

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
gCL = gCL_OpFb_tau1(DelNormArray_OpFb_tau1 == 0.1);
rBall = gCL*dL2norm;
th1Ball4 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball4 = cBall(2) + rBall*sin(aBall);

figure;clf;
patch(th1Ball4,th2Ball4,color1);
hold on;patch(th1Ball1,th2Ball1,color4);
plot(etabar(1)+2*pi,etabar(2),'k',eta0(1)+2*pi, eta0(2),'ko',...
    etaf(1)+2*pi, etaf(2),'kx','LineWidth',2,'MarkerSize',6);axis equal;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;box on;grid on;
set(gca,'xlim',[-0.0612    5.0826]);
set(gca,'ylim',[-3.5532    0.5038]);
alpha(alphaSetting);

%% Plot 5 evrything in same plot

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
gCL = gCL_OpFb_tau1(DelNormArray_OpFb_tau1 == 0.1);
rBall = gCL*dL2norm;
th1Ball4 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball4 = cBall(2) + rBall*sin(aBall);

figure;clf;
patch(th1Ball4,th2Ball4,color1);
hold on;patch(th1Ball3,th2Ball3,color2);
patch(th1Ball2,th2Ball2,color3);
patch(th1Ball1,th2Ball1,color4);
plot(etabar(1)+2*pi,etabar(2),'k',eta0(1)+2*pi, eta0(2),'ko',...
    etaf(1)+2*pi, etaf(2),'kx','LineWidth',2,'MarkerSize',6);axis equal;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;box on;grid on;
set(gca,'xlim',[-0.0612    5.0826]);
set(gca,'ylim',[-3.5532    0.5038]);
alpha(alphaSetting);

%% Plot 6

figure;clf;
tau = 0.1;
[mag,~,wout] = bode(tf([1 0],[tau 1]),{1,1e4});
semilogx(wout,db(mag(:)),color1,'LineWidth',2)

hold on;
tau = 0.05;
[mag,~,wout] = bode(tf([1 0],[tau 1]),{1,1e4});
semilogx(wout,db(mag(:)),color2,'LineWidth',2)

hold on;
tau = 0.02;
[mag,~,wout] = bode(tf([1 0],[tau 1]),{1,1e4});
semilogx(wout,db(mag(:)),color3,'LineWidth',2)

hold on;
tau = 0;
[mag,~,wout] = bode(tf([1 0],[tau 1]),{1,1e4});
semilogx(wout,db(mag(:)),'Color',color4,'LineWidth',2)

xlim([1 1e3])
legend('\tau = 0.1','\tau = 0.05','\tau = 0.02','\tau = 0','Location','northwest','FontSize',14)
title('Bode Diagram');
xlabel('Frequency (rad/s)','FontSize',14);
ylabel('Magnitude (dB)','FontSize',14);
hold off;box on;grid on;
alpha(alphaSetting);