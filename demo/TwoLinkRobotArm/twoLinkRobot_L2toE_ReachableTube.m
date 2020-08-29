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
g = sqrt(max(eig(Woc)));
tTotal = toc(t0);

% Save
save(mfilename);
return;

%% Plot Induced Gains
figure
tvplot(g,'LineWidth',2);
ylabel('Induced L2-to-Euclidean Gain');
box on; grid on;set(gca,'GridLineStyle','--');

%% Plot2D
% Scale to have disturbance of size 5
dNorm = 5;
gScl = g*dNorm;
etabar.InterpolationMethod = 'Linear';
gScl = evalt(gScl,etabar.Time);
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

Tall = linspace(0,Tf,100);
NT = length(Tall);
NBall = 50;
for i = 1:NT
    figure(f3);
    rBall = tvsubs(gScl,Tall(i));
    cBall = tvsubs( etabar(1:2), Tall(i));    
    aBall = linspace(0,2*pi,NBall);
    th1Ball = cBall(1) + rBall*cos(aBall) + 2*pi;
    th2Ball = cBall(2) + rBall*sin(aBall);
    p1 = patch(th1Ball,th2Ball,'c');
    % p1.EdgeColor = [0.5 0.5 0.5];
    
    figure(f4);
    x2Ball = L1*cos(th1Ball) + L2*cos( th1Ball + th2Ball );
    y2Ball = L1*sin(th1Ball) + L2*sin( th1Ball + th2Ball );
    p2 = patch(x2Ball,y2Ball,'c');    
    % p2.EdgeColor = [0.5 0.5 0.5];
end
figure(f3);
plot(theta1f,theta2f,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);
clear alpha, alpha(0.05);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
zlabel('Time (sec)','FontSize',14);
title('$\|d\|_{2,[0,T]} = 5$','Interpreter','latex','FontSize',14);
axis equal;xlim([0 5.5]);ylim([-3.5 0])

figure(f4);
plot(x2f, y2f,'ko','MarkerFaceColor','w','MarkerSize',4);
alpha(0.05);
title('$\|d\|_{2,[0,T]} = 5$','Interpreter','latex','FontSize',14);
xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
zlabel('Time (sec)','FontSize',14);

%% Plot Trajectory: Cartesian Coordinates and Snapshots of Robot Position

figure; hold on;box on;grid on;set(gca,'GridLineStyle','--');
plot(x2,y2,'k','LineWidth', 2.5);
plot(x2f, y2f,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);
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

text(0.01,-0.25,'t=1.5','FontSize',12);
text(-0.4,0.1,'t=2.5','FontSize',12);
text(-0.03,0.3,'t=4.0','FontSize',12);

%% Plot3D
% Scale to have disturbance of size 5
dNorm = 5;
gScl = g*dNorm;
etabar.InterpolationMethod = 'Linear';
gScl = evalt(gScl,etabar.Time);
theta = evalt( etabar(1:2), etabar.Time);
time = etabar.Time;
theta1 = reshapedata(theta(1)+2*pi);
theta2 = reshapedata(theta(2));
theta1f = theta1(end);
theta2f = theta2(end);

f1 = figure;    
hold on;box on;grid on;set(gca,'GridLineStyle','--');
plot3(theta1, theta2, time, 'k','LineWidth', 2.5);

f2 = figure;
hold on;box on;grid on;set(gca,'GridLineStyle','--');
x2 = L1.*cos(theta1) + L2.*cos(theta1+theta2);
y2 = L1.*sin(theta1) + L2.*sin(theta1+theta2);
plot3(x2,y2,time,'k','LineWidth', 2.5);
x2f = L1*cos(theta1f) + L2*cos(theta1f+theta2f);
y2f = L1*sin(theta1f) + L2*sin(theta1f+theta2f);

Tall = linspace(0,Tf,150);
NT = length(Tall);
NBall = 50;
for i = 1:NT
    figure(f1);
    rBall = tvsubs(gScl,Tall(i));
    cBall = tvsubs( etabar(1:2), Tall(i));    
    aBall = linspace(0,2*pi,NBall);
    th1Ball = cBall(1) + rBall*cos(aBall) + 2*pi;
    th2Ball = cBall(2) + rBall*sin(aBall);
    p1 = patch(th1Ball,th2Ball,Tall(i)*ones(NBall,1),'c');
    p1.EdgeColor = [0.5 0.5 0.5];
    
    figure(f2);
    x2Ball = L1*cos(th1Ball) + L2*cos( th1Ball + th2Ball );
    y2Ball = L1*sin(th1Ball) + L2*sin( th1Ball + th2Ball );
    p2 = patch(x2Ball,y2Ball,Tall(i)*ones(NBall,1),'c');    
    p2.EdgeColor = [0.5 0.5 0.5];
end
figure(f1);
plot3(theta1f,theta2f,Tf,'ko','MarkerFaceColor','w','MarkerSize',4,'LineWidth', 1.5);
clear alpha, alpha(0.1);
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
zlabel('Time (sec)','FontSize',14);
title('$\|d\|_{2,[0,T]} = 5$','Interpreter','latex','FontSize',14);

figure(f2);
plot3(x2f, y2f,Tf,'ko','MarkerFaceColor','w','MarkerSize',4);
clear alpha, alpha(0.1);
title('$\|d\|_{2,[0,T]} = 5$','Interpreter','latex','FontSize',14);
xlabel('x (m)','FontSize',14);
ylabel('y (m)','FontSize',14);
zlabel('Time (sec)','FontSize',14);