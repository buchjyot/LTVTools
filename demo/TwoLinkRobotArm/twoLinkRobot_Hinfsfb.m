%% twoLinkRobot_Hinfsfb
% Design a Hinfinity State Feedback Controller & perform analysis

%% Load twoLinkRobot_BuildLTVModel.mat file
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_SpecifyOptions.mat');

%% Finite-Horizon Full-Information Hinf State-Feedback Design
% Full information synthesis
[KSfb,CLsfb,gsfb] = tvhinfsfb(Gsfb,Nu,NE,tvhopt);
fprintf(' Synthesis TVHINFSFB Bisection Bound:');
disp(gsfb);

% Verify bisection results
[gtvn,dWcg] = tvnorm(CLsfb,NE,tvnopt);
fprintf(' TVNORM Bisection Bound:');
disp(gtvn(2));

% Unweight the plant
CLsfb = blkdiag(1/Wu,1/WE)*CLsfb*blkdiag(1/Wd);

%% Nominal Closed-Loop Analysis
fprintf('=================================================\n');
fprintf(' Nominal Closed-Loop Analysis [HINFSFB]\n')
fprintf('=================================================\n');

% Compute L2toE norm for th1 and th2 from d1 and d2
[gtvnE,dWc] = tvnorm(CLsfb(3:4,:),2,tvnopt);
fprintf(' Euclidean TVNORM Closed Loop:');
disp(gtvnE(2));

%% Closed-Loop Linear and Nonlinear Simulations

% Normalized Disturbance
dNorm = 0.1;
dCL = dNorm*dWc/tvnorm(dWc);

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin] = tvlsim(CLsfb(3:4,:),dCL,tvopt);
yf = tvsubs(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dCL);
fprintf(' Nominal CL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin,Tgrid);

% Nonlinear Simulation with Negative Feedback Convention
d = dCL;
model = 'TwoLinkRobotCL_Sfb';
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
title(sprintf('Simulations with ||dCL||=%.1f',dNorm));

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
title(sprintf('Simulations with ||dCL||=%.1f',dNorm))