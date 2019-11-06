%% twoLinkRobot_LQR
% This file designs finite-horizon LQR and analyzes the nominal and robust
% performance

%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat');

%% LQR Design
% Finite-Horizon LQR State-Feedback Design
Q = 100*diag([1 1 0.1 0.1]);
R = 0.1*diag([1,1]);
F = diag([1 1 0.1 0.1]);
Klqr = tvlqr(A,B,Q,R,[],F,[T0,Tf],tvopt);

% Evaluate Klqr and G on the same time grid
[G,Klqr] = evalt(G,Klqr,union(G.Time,Klqr.Time));

% Build closed-loop from disturbance to linearized state x
Tlqr = feedback(G,Klqr);

%% Nominal Closed-Loop Analysis
fprintf('=================================================\n');
fprintf(' Nominal Closed-Loop Analysis [LQR]\n')
fprintf('=================================================\n');

% Nominal Closed-Loop Loop L2 to E Gain
% Expeced 0.05509 (MS Thesis)
[gCLNom,dCL,infoCL] = tvnorm(Tlqr(1:2,:),2,tvnopt);
fprintf(' Bounds on Nominal CL Gain = [%4.4f, %4.4f] \n',...
    gCLNom(1),gCLNom(2));
dNorm = 5;
dCL = dNorm*dCL/tvnorm(dCL);

%% Closed-Loop Linear and Nonlinear Simulations

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin] = tvlsim(Tlqr(1:2,:),dCL,tvopt);
yf = tvsubs(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dCL);
fprintf(' Nominal CL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin,Tgrid);

% Nonlinear Simulation with Negative Feedback Convention
model = 'TwoLinkRobotCL_Sfb';
d = dCL;
KSfb = -Klqr;
sim(model,[T0 Tf]);
yfNL = x(end,1:2);
gNL = norm(yfNL) / tvnorm(dCL);
fprintf(' Nominal CL Gain NonLinear Sim = %4.4f\n',gNL);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid);

% Plot Results
figure;clf
subplot(311)
tvplot(dCL(1),'b',dCL(2),'r--');
ylabel('d (Nm)');
legend('d1','d2');
title(sprintf('Simulations with ||dCL||=%.3f',dNorm));

subplot(312)
tvplot(etabar(1),'k-.', etaNL(1),'r', etaLin(1),'b--');
ylabel('\theta_1 (rad)');
legend('Trim','NL','Lin');

subplot(313)
tvplot(etabar(2),'k-.', etaNL(2),'r', etaLin(2),'b--');
xlabel('Time (sec)');
ylabel('\theta_2 (rad)');

figure;clf;
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title(sprintf('Simulations with ||dCL||=%.3f',dNorm))

%% Robust Closed-Loop Analysis
fprintf('=================================================\n');
fprintf(' Robust Closed-Loop Analysis [LQR] \n')
fprintf('=================================================\n');

% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
systemnames = 'G Klqr';
inputvar = '[w; d(2)]';
outputvar = '[d(2)-Klqr(2) ; G(1:2)]';
input_to_G = '[d(1)-Klqr(1); d(2)-Klqr(2)+w]';
input_to_Klqr = '[G]';
cleanupsysic = 'yes';
Tnom = sysic;

% Uncertainity Scalling
Lscl = diag([sqrt(DelNorm) 1 1]);
Rscl = diag([sqrt(DelNorm) 1 1]);
Tnom = Lscl*Tnom*Rscl;

% Robust (worst-case) L2 to E Gain
% Expected 0.062 (MS Thesis)
% Deprecated Syntax: [gCL,wcinfoCL] = tvrobL2toE(Tnom,v,p,Tf,tlmi,tSp);
Tnom.UserData = IQCParam;
[gCL,wcinfoCL] = tvwcgain(Tnom,2,tvwcopt);

%% Evaluate closed-loop with a specific bad perturbation

% Wrap in "Bad" Perturbation
% The *.Mat file contains 100 random samples of perturbations (with
% norm ||Delta||=0.8). The code below selects a perturbation and
% normalizes it to have ||Delta||=1.
% [Note: The selected Delta is the one that maximizes the closed-loop gain
%  out of these 100 random samples.  The calculation of the closed-loop
%  norms on these 100 samples is skipped to reduce computation.]
load('RandomDelta.mat');
DeltaBad = DeltaSim{11} / norm(DeltaSim{11},inf);
Tbad = lft( DeltaBad ,Tnom);

% Evaluate closed-loop gain of Tbad
% Gain lower bound should be close to 0.0575
fprintf('\n ---- Evaluate "Bad" Perturbation \n');
[gWC,dWC] = tvnorm(Tbad,2,tvnopt);
fprintf('\n Closed-loop gain with worst-case Delta = %4.4f',gWC(1))

% Simulate linear system and evaluate gain
yBad = tvlsim(Tbad,dWC,tvopt);
gWC2 = norm( tvsubs(yBad,Tf) ) / tvnorm(dWC);
fprintf('\n Closed-loop gain with worst-case Delta and dist. = %4.4f\n',gWC2)

%% Simulate and Plot Trajectories / Norm Bound

% Disturbances are scaled to have norm ||d|| = dL2norm
dL2norm = 5;

% e(T) ball: Center at final (trim) angles with radius gCL*dL2norm
rBall = gCL*dL2norm;
cBall = tvsubs( etabar(1:2), Tf);
NBall = 50;
aBall = linspace(0,2*pi,NBall);
th1Ball = cBall(1) + rBall*cos(aBall) + 2*pi;
th2Ball = cBall(2) + rBall*sin(aBall);

% Convert e(T) ball to Cartesian coordinates
x2Ball = L1*cos(th1Ball) + L2*cos( th1Ball + th2Ball );
y2Ball = L1*sin(th1Ball) + L2*sin( th1Ball + th2Ball );

% Draw e(T) ball and corresponding Cartesian coordinate region
f1 = figure;
clf;
axis([0 5, -3.5 0.5]);
patch(th1Ball,th2Ball,'c');
hold on;

f2 = figure;
clf;
axiswidth = L1 + L2 + 0.2*(L1 + L2);
axis equal;
axis(axiswidth*[-1, 1, -1 ,1]);hold on;
patch(x2Ball,y2Ball,'c');
hold on;

% Load Random Disturbances and append worst-case disturbance
load('RandomDisturbances.mat');
d = [d, {dWC}];

% Simulate linearized system
for i = 1:numel(d)
    % Display count i
    if floor(i/10)==ceil(i/10)
        fprintf('\n i=%d ',i)
    end
    
    % Disturbances are scaled to have norm ||d||=dL2norm
    di = d{i};
    di = di*dL2norm/tvnorm(di);
    
    % Simulate linear system
    yi = tvlsim(Tbad,di,tvopt);
    
    % Superimpose the trim trajectory to obtain actual angles
    % (See twoLinkRobot_BuildLTVModel for trim trajectory construction)
    theta = evalt( etabar(1:2), yi.Time);
    theta.InterpolationMethod = 'Linear';
    theta = yi+ theta;
    theta1 = theta(1)+2*pi;
    theta2 = theta(2);
    
    % Plot trajectory and mark final point
    figure(f1);
    plot(theta1,theta2);
    theta1f = tvsubs(theta1,Tf);
    theta2f = tvsubs(theta2,Tf);
    plot(theta1f,theta2f,'ko','MarkerFaceColor','w');
    
    % Cartesian coordinates of link 2 tip
    x2 = L1*cos(theta1) + L2*cos(theta1+theta2);
    y2 = L1*sin(theta1) + L2*sin(theta1+theta2);
    
    % Plot trajectory and mark final point
    figure(f2);
    plot(x2,y2);
    x2f = tvsubs(x2,Tf);
    y2f = tvsubs(y2,Tf);
    plot(x2f,y2f,'ko','MarkerFaceColor','w');
end
fprintf('\n');

% Add labels and legends
figure(f1);
grid on;
xlabel('\theta_1 (rads)');
ylabel('\theta_2 (rads)');
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;

figure(f2);
grid on;
xlabel('x (m)'); ylabel('y (m)');
title(['Closed-Loop with ||d|| <= ' num2str(dL2norm)]);
hold off;