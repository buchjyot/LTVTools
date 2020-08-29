%% twoLinkRobot_ReachableTube.m

% Load Data
load('twoLinkRobot_LQR.mat','Tlqr');
load('twoLinkRobot_BuildLTVModel.mat','T0','Tf','L1','L2','etabar');

% Get LTV StateSpace Matrices
SYS = Tlqr(1:2,:);
[A,B,C] = ssdata(SYS);

%% Uncertainty in Initial Condition in theta1, theta2 space
E0_CASE = 2;
switch E0_CASE
    case 1
        rng(5757);
        M0 = randn(2);
        E0 = M0'*M0;
        E0 = E0/norm(E0);
        E0 = E0*rad2deg(1);
    case 2
        E0 = 5*eye(2);
end

%% Plot Trajectory
f1 = figure;clf;
theta = evalt( etabar(1:2), etabar.Time);
time = etabar.Time;
theta1 = theta(1)+2*pi;
theta2 = theta(2);
theta10 = tvsubs(theta1,T0);
theta20 = tvsubs(theta2,T0);
theta1f = tvsubs(theta1,Tf);
theta2f = tvsubs(theta2,Tf);
hold on;box on;grid on;set(gca,'GridLineStyle','--');
plot(theta1, theta2, 'k','LineWidth', 2.5);
plot(theta1f,theta2f,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
axis equal;xlim([0 5.5]);ylim([-3.5 0])
p1 = plot_gaussian_ellipsoid([theta10,theta20],inv(E0));
patch(p1.XData,p1.YData,'b');
alpha(0.1);

%% Start Power-Iterations
E04d = blkdiag(E0,diag(inf(2,1)));
pSpec = poweritSignalSpec('NE',2,'InitialConditions','free','InitialCondCostMat',E04d,'InputL2Norm',5);
pOpt = poweritOptions('Display','on','StoreAllIter',true);
[glb,dwc,info] = powerit(SYS,[T0,Tf],pSpec,pOpt);

%% Plot2D
% Scale to have disturbance of size 5
etabar.InterpolationMethod = 'Linear';

% Copy first figure
figure(f1);
ax = gca;
f2 = figure;
a2 = copyobj(ax,f2);

% Worst-case trajectory
Ywc = info.Ywc;
etabar = evalt(etabar,Ywc.Time);
theta1wc = etabar(1) + Ywc(1) + 2*pi;
theta2wc = etabar(2) + Ywc(2);
plot(theta1wc, theta2wc, 'r-.','LineWidth', 2.5);
theta1fwc = tvsubs(theta1wc,Tf);
theta2fwc = tvsubs(theta2wc,Tf);
plot(theta1fwc,theta2fwc,'ro','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);

% Plot the ball
rBall = norm([theta1fwc;theta2fwc] - [theta1f;theta2f]);
[XData,YData] = getCircleData([theta1f,theta2f],rBall);
patch(XData,YData,'r');
alpha(0.1);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);

%% Monte-Carlo Sims
figure(f2);
ax = gca;
f3 = figure;
a3 = copyobj(ax,f3);hold on;

load('RandomDisturbances.mat','d');
rng(0);
Nd = length(d);
Y = cell(Nd,1);
x0a = cell(Nd,1);
theta1i = cell(Nd,1);
theta1fi = cell(Nd,1);
theta2fi = cell(Nd,1);
theta2i = cell(Nd,1);
for i = 1:Nd
    % Linear Sim
    x0 = randn(2,1);
    nx0 = x0'*E0*x0;
    x0a{i} = (rand*x0)/sqrt(nx0);
    Y{i} = tvlsim(SYS,d{i},[T0,Tf],[x0a{i};zeros(2,1)]);    
    
    % Trajectory
    etabari = evalt(etabar,Y{i}.Time);
    theta1i{i} = etabari(1) + Y{i}(1) + 2*pi;
    theta2i{i} = etabari(2) + Y{i}(2);
    plot(theta1i{i}, theta2i{i},'b');
    
    % Final Position
    theta1fi{i} = tvsubs(theta1i{i},Tf);
    theta2fi{i} = tvsubs(theta2i{i},Tf);
    plot(theta1fi{i},theta2fi{i},'bo','MarkerFaceColor','w','MarkerSize',4);
end