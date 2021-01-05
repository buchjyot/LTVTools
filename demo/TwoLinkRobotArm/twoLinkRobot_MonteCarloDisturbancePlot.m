%% twoLinkRobot_LQR
% This file designs finite-horizon LQR and analyzes the nominal and robust
% performance

%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_SpecifyOptions.mat');

%% LQR Design
% Finite-Horizon LQR State-Feedback Design
Q = 100*diag([1 1 0.1 0.1]);
R = 0.1*diag([1,1]);
F = diag([1 1 0.1 0.1]);
Klqr = tvlqr(A,B,Q,R,[],F,[T0,Tf],tvopt);

%% Nonlinear Monte-Carlo Sims with Sampled Disturbances

% Model
model = 'TwoLinkRobotCL_Sfb';
KSfb = -Klqr;
load_system(model);

% Disturbances are scaled to have norm ||d|| = dL2norm
dNorm = 5;
f1 = figure;
clf;
axis([0 5, -3.5 0.5]);
hold on;grid on;box on;

f2 = figure;
clf;
axiswidth = L1 + L2 + 0.2*(L1 + L2);
axis equal;
axis([-0.6 0.6 -0.3 0.5]);hold on;
hold on;grid on; box on;

% Load Random Disturbances and append worst-case disturbance
load('RandomDisturbances.mat');
dall = [d {tvmat(zeros(2,1),[T0:1:Tf])}];
Nd = numel(dall);
tvopt1 = tvopt;
tvopt1.OdeSolver = 'ode45';

% Simulate linearized system
theta1f = zeros(Nd,1);
theta2f = zeros(Nd,1);
x2f = zeros(Nd,1);
y2f = zeros(Nd,1);
for i = 1:Nd
    % Display count i
    if floor(i/10)==ceil(i/10)
        fprintf('\n i=%d ',i)
    end
    
    % Disturbances are scaled to have norm ||d||=dL2norm
    di = dall{i};
    if i~=Nd
        di = di*dNorm/tvnorm(di);
    end
    
    % Nonlinear Simulation with Negative Feedback Convention
    d = di;
    sim(model,[T0 Tf]);
    
    % Superimpose the trim trajectory to obtain actual angles
    % (See twoLinkRobot_BuildLTVModel for trim trajectory construction)
    theta  = tvmat(eta(:,1:2)',tsim);
    theta1 = theta(1)+2*pi;
    theta2 = theta(2);
    
    % Plot trajectory and mark final point
    figure(f1);
    if i == Nd
        plot(theta1,theta2,'k','LineWidth',2.5);
    else
        plot(theta1,theta2,'c','LineWidth',1);
    end
    theta1f(i) = tvsubs(theta1,Tf);
    theta2f(i) = tvsubs(theta2,Tf);    
    
    % Cartesian coordinates of link 2 tip
    x2 = L1*cos(theta1) + L2*cos(theta1+theta2);
    y2 = L1*sin(theta1) + L2*sin(theta1+theta2);
    
    % Plot trajectory and mark final point
    figure(f2);
    if i == Nd
        plot(x2,y2,'k','LineWidth',2.5);
    else
        plot(x2,y2,'c','LineWidth',1);
    end
    x2f(i) = tvsubs(x2,Tf);
    y2f(i) = tvsubs(y2,Tf);    
end

% Plot end markers
for i = 1:Nd
    figure(f1);hold on;
    plot(theta1f(i),theta2f(i),'ko','MarkerFaceColor','w','MarkerSize',5);
    
    figure(f2);hold on;
    plot(x2f(i),y2f(i),'ko','MarkerFaceColor','w','MarkerSize',5);
end

% Add labels and legends
figure(f1);
grid on;
xlabel('\theta_1 (rads)','FontSize',14);
ylabel('\theta_2 (rads)','FontSize',14);
hold off;

figure(f2);
grid on;
xlabel('x (m)','FontSize',14); ylabel('y (m)','FontSize',14);
hold off;