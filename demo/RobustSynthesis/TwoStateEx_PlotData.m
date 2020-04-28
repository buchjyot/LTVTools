%% TwoStateEx_UncSweep
load('TwoStateEx_UncSweep.mat');
load('TwoStateEx.mat','beta');
figure;clf;grid on;box on;hold on;
plot(UL,wcgain1,'-ob',UL,wcgain2,'-.rs','LineWidth',2);
legend('$\tilde{T}_0$',['$\tilde{T}_{' num2str(beta) '}$'],...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)');
ylabel('Worst-Case Gain');
XL = xlim;
YL = ylim;

figure;clf;grid on;box on;hold on;
plot(UL(1),wcgain1(1),'-ob',UL(7),wcgain2(7),'-.rs','LineWidth',2);
legend('$\tilde{T}_0$',['$\tilde{T}_{' num2str(beta) '}$'],...
    'interpreter','latex','location','northwest','fontsize',14);
xlabel('Uncertainty Level (\beta)');
ylabel('Worst-Case Gain');
xlim(XL);ylim(YL);
clear;

%% TwoStateEx_L2toE_NominalAnalysis
load('TwoStateEx_L2toE_NominalAnalysis.mat')
figure;clf;
X0n = tvsubs(XCLn,T0);
X0r = tvsubs(XCLr,T0);
XTfn = tvsubs(XCLn,Tf);
XTfr = tvsubs(XCLr,Tf);
grid on;hold on;box on;
[X1,Y1] = getCircleData(X0n',gCLn(1)*dScl);
[X2,Y2] = getCircleData(X0r',gCLr(1)*dScl);
plot(X1,Y1,'b','LineWidth',2);
plot(X2,Y2,'r','LineWidth',2);
tvplot(XCLn(1),XCLn(2),'-.b',XTfn(1),XTfn(2),'bo','MarkerFaceColor','b','MarkerSize',6,'LineWidth',2);
tvplot(XCLr(1),XCLr(2),'-.r',XTfr(1),XTfr(2),'ro','MarkerFaceColor','r','MarkerSize',6,'LineWidth',2);
plot(X0n(1),X0n(2),'kx','MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
legend('$T_{0 (d_1\rightarrow e_E)}$',...
    ['$T_{' num2str(beta) '(d_1\rightarrow e_E)}$'],'interpreter','latex','orientation','horizontal','fontsize',14,...
    'location','northoutside')
axis equal
yl = ylim;
ylim(yl+0.1*yl);
xlim([-0.6 0.6])
clear;

%% TwoStateEx_L2toE_RobustAnalysis
load('TwoStateEx_L2toE_RobustAnalysis.mat');
figure;clf;
X0n = tvsubs(XCLn,T0);
X0r = tvsubs(XCLr,T0);
XTfn = tvsubs(XCLn,Tf);
XTfr = tvsubs(XCLr,Tf);
grid on;hold on;box on;
[X1,Y1] = getCircleData([0 0],wcgUB1*dScl);
[X2,Y2] = getCircleData([0 0],wcgUB2*dScl);
plot(X1,Y1,'b','LineWidth',2);
plot(X2,Y2,'r','LineWidth',2);
p1 = tvplot(XCLn(1),XCLn(2),'-.b',XTfn(1),XTfn(2),'bo','MarkerFaceColor','b','MarkerSize',6,'LineWidth',2);
p2 = tvplot(XCLr(1),XCLr(2),'-.r',XTfr(1),XTfr(2),'ro','MarkerFaceColor','r','MarkerSize',6,'LineWidth',2);
plot(X0n(1),X0n(2),'kx','MarkerEdgeColor','k','MarkerFaceColor','k');
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
legend('$T_{0 (d_1\rightarrow e_E)}$','$T_{0.6 (d_1\rightarrow e_E)}$','interpreter','latex','fontsize',14,...
    'orientation','horizontal','location','northoutside')
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
axis equal
yl = ylim;
ylim(yl+0.1*yl)
xlim([-1 1])
clear

%% Plot Iteration Progress
load('TwoStateEx_RobustSynthesis.mat','robinfo');
figure;clf;
wcg = cellfun(@(x) x.wcgain,robinfo);
plot(2:20,wcg(2:20),'k-*','LineWidth',1.5,'MarkerSize',8);
ylabel('Worst-Case Gain');
xlabel('Iteration Count (i)');
grid on;box on;
clear