%% twoLinkRobot_OpenLoop
% OpenLoop Analysis for Linearized Two-Link Robot Model

%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_SpecifyOptions.mat');

%% Nominal Open-Loop Analysis
fprintf('==================================================\n');
fprintf('Nominal Open-Loop Analysis\n')
fprintf('==================================================\n');

% Nominal Open Loop L2 to E Gain
% Expected 528.9 (MS Thesis)
NE = 2;
tvnopt.Bounds = [0 1e3];
[gOLNom,dOL,infoOL] = tvnorm(G(1:2,:),NE,tvnopt);
fprintf(' Bounds on Nominal OL Gain = [%4.4f, %4.4f] \n',...
    gOLNom(1),gOLNom(2));
dNorm = 0.001;
dOL = dNorm*dOL/tvnorm(dOL);

%% Open-Loop Linear and Nonlinear Simulations

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin] = tvlsim(G(1:2,:),dOL,tvnopt);
yf = tvsubs(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dOL);
fprintf(' Nominal OL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin,Tgrid);

% Nonlinear Simulation
d = dOL;
model = 'TwoLinkRobotOL';
sim(model,[T0 Tf]);
yfOLNL = x(end,1:2);
gNL = norm(yfOLNL) / tvnorm(dOL);
fprintf(' Nominal OL Gain NonLinear Sim = %4.4f\n',gNL);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid);

% Plot Results
figure;clf;
subplot(311)
tvplot(dOL(1),'b',dOL(2),'r--');
ylabel('d (Nm)');
legend('d1','d2');
title(sprintf('Simulations with ||dOL||=%.1f',dNorm));

subplot(312)
tvplot(etabar(1),'k-.',etaNL(1),'r',etaLin(1),'b--');
ylabel('\theta_1 (rad)');
legend('Trim','NL','Lin');

subplot(313)
tvplot(etabar(2),'k-.',etaNL(2),'r',etaLin(2),'b--');
xlabel('Time (sec)');
ylabel('\theta_2 (rad)');

figure;clf;
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title(sprintf('Simulations with ||dOL||=%.4f',dNorm))

%% Robust Open-Loop Analysis
fprintf('==================================================\n');
fprintf('Robust Open-Loop Analysis\n')
fprintf('==================================================\n');
tvwcopt.RDEOptions.Bounds = [0 1e3];
if true
    % Uncertain OL System:  Gunc = Fu(Gnom,Delta)
    % Delta is unit norm-bounded LTI uncertainty.  Gnom is constructed
    % so that uncertainty enters into second input channel as 0.8*Delta.
    Gnom = [0 0 sqrt(DelNorm); G(1:2,:)*[0 1 0; sqrt(DelNorm) 0 1] ];
    
    % Robust (worst-case) L2 to E Gain
    % Expected 995 (MS Thesis)
    % Deprecated Syntax: [gOL,wcinfoOL] = tvrobL2toE(Gnom,v,p,Tf,tlmi,tSp);
    [gOL,wcinfoOL] = tvwcgain(Gnom,Delta,NE,tvwcopt);
end