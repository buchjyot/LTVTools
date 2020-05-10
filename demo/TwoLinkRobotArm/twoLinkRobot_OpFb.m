%% twoLinkRobot_OpFb
% Output feedback design using LQR gains and highpass filter to obtain rate
% estimates.

%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_SpecifyOptions.mat');

%% Output-Feedback with Rate Estimates

% High-Pass Filter Design
% High-Pass Filter [Dirty-Derivative] at the output
%
%   Hpf(s) = s/(tau*s+1)
%
% The goal is to estimate angular rates using (possibly noisy) angular
% position measurements

% Finite-Horizon LQR State-Feedback Design
Q = 100*diag([1 1 0.1 0.1]);
R = 0.1*diag([1,1]);
F = diag([1 1 0.1 0.1]);
Klqr = tvlqr(A,B,Q,R,[],F,[T0,Tf],tvopt);

% Evaluate Klqr and G on the same time grid
[G,Klqr] = evalt(G,Klqr,union(G.Time,Klqr.Time));

% Time-Constant
tauHpf = 0.1;

% State Space Realization is used here for above transfer function
Hpf = ss(-1/tauHpf,-1/tauHpf,1/tauHpf,1/tauHpf);

% First Get the Linearized plant with output being th1 and th2
GOpFb = G(1:2,:);

% Output Feedback Controller We have states as [th1 th2 th1dot th2dot]'
KOpFb = Klqr*[1 0;0 1;Hpf 0;0 Hpf];

% Build closed-loop from disturbance to linearized state x
% Total 6 states = 4 Plant + 2 Controller
[GOpFb,KOpFb] = evalt(GOpFb,KOpFb,union(KOpFb.Time,GOpFb.Time));
TOpFb = feedback(GOpFb,KOpFb);

%% Nominal Closed-Loop Analysis With Output Feedback
fprintf('==================================================\n');
fprintf('Nominal Closed-Loop Analysis [Output Feedback]\n')
fprintf('==================================================\n');

% Nominal Closed-Loop Loop L2 to E Gain
[gCLNom,dCL,infoCLOpFb] = tvnorm(TOpFb,2,tvnopt);
fprintf(' Bounds on Nominal CL Gain = [%4.4f, %4.4f] \n',...
    gCLNom(1),gCLNom(2));
dNorm = 5;
dCL = dNorm*dCL/tvnorm(dCL);

%% Closed-Loop Linear and Nonlinear Simulations

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin] = tvlsim(TOpFb,dCL,tvspt);
yf = tvsubs(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dCL);
fprintf(' Nominal CL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin(1:4,:),Tgrid);

% Nonlinear Simulation
d = dCL;
KSfb = -Klqr;
model = 'TwoLinkRobotCL_OpFbHpf';
sim(model,[T0 Tf]);
yfNL = y(end,:);
gNL = norm(yfNL) / tvnorm(dCL);
fprintf(' Nominal CL Gain NonLinear Sim = %4.4f\n',gNL);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid);

% Plot Results
figure;clf;
subplot(311);
tvplot(dCL(1),'b',dCL(2),'r--');
ylabel('d (Nm)');
legend('d1','d2');
title('Simulations with ||dCL||=20')

subplot(312);
tvplot(etabar(1),'k-.',etaNL(1),'r',etaLin(1),'b--');
ylabel('\theta_1 (rad)');
legend('Trim','NL','Lin');

subplot(313);
tvplot(etabar(2),'k-.',etaNL(2),'r',etaLin(2),'b--');
xlabel('Time (sec)');
ylabel('\theta_2 (rad)');

figure;clf;
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title(sprintf('Simulations with ||dCL||=%.3f',dNorm))

%% Robust Closed-Loop Analysis [Output Feedback]
fprintf('==================================================\n');
fprintf('Robust Closed-Loop Analysis [Output Feedback]\n')
fprintf('==================================================\n');

% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
systemnames = 'GOpFb KOpFb';
inputvar = '[w; d(2)]';
outputvar = '[d(2)-KOpFb(2); GOpFb]';
input_to_GOpFb = '[d(1)-KOpFb(1); d(2)-KOpFb(2)+w]';
input_to_KOpFb = '[GOpFb]';
cleanupsysic = 'yes';
Tnom = sysic;

% Uncertainty Norm Bound
DelNorm = 0.1;

% Uncertainity Scalling
Lscl = diag([sqrt(DelNorm) 1 1]);
Rscl = diag([sqrt(DelNorm) 1 1]);
Tnom = Lscl*Tnom*Rscl;

% Robust (worst-case) L2 to E Gain
% Deprecated Syntax: [gCL,wcinfoCL] = tvrobL2toE(Tnom,v,p,Tf,tlmi,tSp);
[gCL,wcinfoCL] = tvwcgain(Tnom,Delta,2,tvwcopt);