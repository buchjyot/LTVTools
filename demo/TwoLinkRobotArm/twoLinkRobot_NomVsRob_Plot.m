%% Load Data
% This file plots the nominal trajectory and comparision for nominal and
% robust euclidean gains
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_SpecifyOptions.mat');
load('twoLinkRobot_HinfDesign.mat','dScl','DelNorm');
clear alpha

% Define center as the final position of the tip of the arm
cBall = tvsubs( etabar(1:2), etabar.Time(end));
NBall = 50;
aBall = linspace(0,2*pi,NBall);
alphaPatch = 0.15;

% Colors
DarkGreen = [0 0.3906 0];
DarkRed = [0.5430 0 0];
%% Plot Nominal Trajectory
% Create figure
figure;hold on;

% Plot Nominal Trajectory
plot(2*pi+etabar(1), etabar(2), 'k','LineWidth',2.5);
plot(2*pi+etaf(1), etaf(2),'ko','MarkerSize',8,'MarkerFaceColor','k');

% Add legends
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
grid on;axis equal;box on;
text(1,-0.65,'t=5','FontSize',14);
text(4,-1.15,'t=0','FontSize',14);
axis([0.5 5, -3.5 0]);

%% Plot Trajectory: Cartesian Coordinates and Snapshots of Robot Position
figure; clf; hold on;
plot( x2bar, y2bar, 'k', 'LineWidth',2.5);
plot( tvsubs(x2bar,Tf), tvsubs(y2bar,Tf),'ko','MarkerSize',8,'MarkerFaceColor','k')
xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
grid on;

hold on;
for ti = [1.5 2.5 4]
    xi = tvsubs([x1bar;x2bar;y1bar; y2bar],ti);
    plot([0; xi(1:2)],[0; xi(3:4)],'g','LineWidth',3);
    %plot(xi(1),xi(3),'ko','MarkerSize',4,'MarkerFaceColor','g');
    plot(xi(2),xi(4),'ko','MarkerSize',8,'MarkerFaceColor','m');
    plot(0,0,'bo','MarkerSize',8,'MarkerFaceColor',DarkGreen);
end
hold off;
axis equal;xlim([-0.6 0.6])

text(0.01,-0.3,'t=1.5','FontSize',12);
text(-0.4,0.1,'t=2.5','FontSize',12);
text(-0.03,0.3,'t=4.0','FontSize',12);
box on;

%% Nominal Synthesis

% Load Data
load('twoLinkRobot_L2toE_NominalAnalysis.mat');

% Create figure
figure;hold on;

% Plot the patch of Nominal Controller performance upper bound
rBall = dScl*gCLn(2);
th1Ball1 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball1 = cBall(2) + rBall*sin(aBall);
p1 = patch(th1Ball1,th2Ball1,'b');
alpha(p1,alphaPatch);

% Plot the patch of Robust Controller performance upper bound
rBal2 = dScl*gCLr(2);
th1Ball2 = cBall(1) + rBal2*cos(aBall) + 2*pi;
th2Ball2 = cBall(2) + rBal2*sin(aBall);
p2 = patch(th1Ball2,th2Ball2,'r');
alpha(p2,alphaPatch);

% Plot Nominal Trajectory
plot(2*pi+etabar(1), etabar(2), 'k','LineWidth',2.5);
plot(2*pi+etaf(1), etaf(2),'ko','MarkerSize',8,'MarkerFaceColor','k');

% Plot worst-case trajectory for nominal controller
XCLn = evalt(XCLn,etabar.Time);
XCLn.InterpolationMethod = 'Spline';
XCLnf = tvsubs(XCLn,XCLn.Time(end));
plot(etabar(1)+2*pi+XCLn(1),etabar(2)+XCLn(2),'--b','LineWidth',2.5);
plot(2*pi+etaf(1)+XCLnf(1), etaf(2)+XCLnf(2),'bo','MarkerFaceColor','b','MarkerSize',8);

% Plot worst-case trajectory for robust controller
XCLr = evalt(XCLr,etabar.Time);
XCLr.InterpolationMethod = 'Spline';
XCLrf = tvsubs(XCLr,XCLr.Time(end));
plot(etabar(1)+2*pi+XCLr(1),etabar(2)+XCLr(2),'--r','LineWidth',2.5);
plot(2*pi+etaf(1)+XCLrf(1), etaf(2)+XCLrf(2),'ro','MarkerFaceColor','r','MarkerSize',8);

% Add legends
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
grid on;axis equal;box on;
axis([0 5, -3.5 0.5]);
title(' ')

ax1 = gca;
f1 = figure;
copyobj(ax1,f1);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
axis([0.2993    1.3978   -1.2404   -0.1877]);
legend('$T_{0 (d\rightarrow e_E)}$',...
    ['$T_{' num2str(0.8) '(d\rightarrow e_E)}$'],'interpreter','latex','orientation','horizontal','fontsize',14,...
    'location','north')
title(' ')

%% Robust Synthesis

% Load data
load('twoLinkRobot_L2toE_RobustAnalysis.mat');

% Create figure
figure;hold on;

% Plot the patch of Nominal Controller performance upper bound
rBall = dScl*wcgUB1;
th1Ball1 = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball1 = cBall(2) + rBall*sin(aBall);
p1 = patch(th1Ball1,th2Ball1,'b');
alpha(p1,alphaPatch);

% Plot the patch of Robust Controller performance upper bound
rBal2 = dScl*wcgUB2;
th1Ball2 = cBall(1) + rBal2*cos(aBall) + 2*pi;
th2Ball2 = cBall(2) + rBal2*sin(aBall);
p2 = patch(th1Ball2,th2Ball2,'r');
alpha(p2,alphaPatch);

% Plot Nominal Trajectory
plot(2*pi+etabar(1), etabar(2), 'k','LineWidth',2.5);
plot(2*pi+etaf(1), etaf(2),'ko','MarkerSize',8,'MarkerFaceColor','k');

% Plot worst-case trajectory for nominal controller
XCLn = evalt(XCLn,etabar.Time);
XCLn.InterpolationMethod = 'Spline';
XCLnf = tvsubs(XCLn,XCLn.Time(end));
plot(etabar(1)+2*pi+XCLn(1),etabar(2)+XCLn(2),'--b','LineWidth',2.5);
plot(2*pi+etaf(1)+XCLnf(1), etaf(2)+XCLnf(2),'bo','MarkerFaceColor','b','MarkerSize',8);

% Plot worst-case trajectory for robust controller
XCLr = evalt(XCLr,etabar.Time);
XCLr.InterpolationMethod = 'Spline';
XCLrf = tvsubs(XCLr,XCLr.Time(end));
plot(etabar(1)+2*pi+XCLr(1),etabar(2)+XCLr(2),'--r','LineWidth',2.5);
plot(2*pi+etaf(1)+XCLrf(1), etaf(2)+XCLrf(2),'ro','MarkerFaceColor','r','MarkerSize',8);

% Add legends
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
grid on;axis equal;box on;
axis([0 5, -3.5 0.5]);
title(' ')

ax2 = gca;
f2 = figure;
copyobj(ax2,f2);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
title('')
axis([0.1124    1.4967   -1.3561   -0.0295]);
legend('$T_{0 (d\rightarrow e_E)}$',...
    ['$T_{' num2str(0.8) '(d\rightarrow e_E)}$'],'interpreter','latex','orientation','horizontal','fontsize',14,...
    'location','north')

%% TwoStateEx_UncSweep
load('twoLinkRobot_UncSweep.mat');
figure;clf;grid on;box on;hold on;
plot(UL,wcgain1,'-ob',UL,wcgain2,'-.rs','LineWidth',2);
legend('$\tilde{T}_0$',['$\tilde{T}_{' num2str(DelNorm) '}$'],...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
XL = xlim;
YL = ylim;

figure;clf;grid on;box on;hold on;
plot(UL(1),wcgain1(1),'-ob',UL(9),wcgain2(9),'-.rs','LineWidth',2);
legend('$\tilde{T}_0$',['$\tilde{T}_{' num2str(DelNorm) '}$'],...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)','FontSize',14);
ylabel('Worst-Case Gain','FontSize',14);
xlim(XL);ylim(YL);
clear;