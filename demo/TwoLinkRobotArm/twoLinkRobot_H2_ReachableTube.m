%% twoLinkRobot_ReachableTube.m

% Load Data
load('twoLinkRobot_LQR.mat','Tlqr');
load('twoLinkRobot_BuildLTVModel.mat','T0','Tf','L1','L2','etabar');

% Get LTV StateSpace Matrices
[A,B,C] = ssdata(Tlqr(1:2,:));

% Compute Gramian and Maximum Eigenvalue
t0 = tic;
Wc = cdle(A,B,zeros(4),C.Time);
Woc = C*Wc*C';
g = trace(Woc);
tTotal = toc(t0);

% Save
save(mfilename);
return;

%% Plot Induced Gains
figure
tvplot(g,'LineWidth',2);
ylabel('H2 Norm');
box on; grid on;set(gca,'GridLineStyle','--');

%% Plot2D
% Scale to have disturbance of size 5
etabar.InterpolationMethod = 'Linear';
g = evalt(g,etabar.Time);
theta = evalt( etabar(1:2), etabar.Time);
time = etabar.Time;
theta1 = theta(1)+2*pi;
theta2 = theta(2);
theta1f = tvsubs(theta1,Tf);
theta2f = tvsubs(theta2,Tf);

f3 = figure;
hold on;box on;grid on;set(gca,'GridLineStyle','--');
plot(theta1, theta2, 'k','LineWidth', 2.5);

f4 = figure;
hold on;box on;grid on;set(gca,'GridLineStyle','--');
x2 = L1.*cos(theta1) + L2.*cos(theta1+theta2);
y2 = L1.*sin(theta1) + L2.*sin(theta1+theta2);
plot(x2,y2,'k','LineWidth', 2.5);
x2f = L1*cos(theta1f) + L2*cos(theta1f+theta2f);
y2f = L1*sin(theta1f) + L2*sin(theta1f+theta2f);

Tall = linspace(0,Tf,150);
NT = length(Tall);
NBall = 50;
for i = 1:NT
    figure(f3);
    cBall = tvsubs( [theta1;theta2], Tall(i));
    WocT = tvsubs(Woc, Tall(i));
    h1 = plot_gaussian_ellipsoid(cBall,WocT,3);
    th1Ball = h1.XData;
    th2Ball = h1.YData;
    p1 = patch(th1Ball,th2Ball,'c');
    h1.Color = 'k';
    
    figure(f4);
    x2Ball = L1*cos(th1Ball) + L2*cos( th1Ball + th2Ball );
    y2Ball = L1*sin(th1Ball) + L2*sin( th1Ball + th2Ball );
    p2 = patch(x2Ball,y2Ball,'c');
    p2.EdgeColor = [0.5 0.5 0.5];
end
figure(f3);
plot(theta1f,theta2f,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);
clear alpha, alpha(0.05);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
zlabel('Time (sec)','FontSize',14);
axis equal;xlim([0 5.5]);ylim([-3.5 0])

figure(f4);
plot(x2f, y2f,'ko','MarkerFaceColor','w','MarkerSize',4);
alpha(0.05);
xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
zlabel('Time (sec)','FontSize',14);