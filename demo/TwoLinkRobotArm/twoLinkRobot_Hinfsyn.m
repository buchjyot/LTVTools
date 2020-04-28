%% twoLinkRobot_Hinfsyn
% Design Hinfinity Output Feedback Controller  & perform analysis

%% Load twoLinkRobot_BuildLTVModel.mat file
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_SpecifyOptions.mat');
tvnopt.RelTol = 5e-3;

%% Design Output Feedback Hinf Controller

% Synthesis Step
[Knom,CLnom,gnom,nominfo] = tvhinfsyn(Gnom,Ny,Nu,NE,tvhopt);
fprintf(' Synthesis TVHINFSYN Bound:');
disp(gnom);

% Verify bisection results
[g,dwc] = tvnorm(CLnom,NE,tvnopt);
fprintf(' TVNORM Bisection Bound:');
disp(g(2));

%% Nominal Closed-Loop Analysis

fprintf('=================================================\n');
fprintf(' Nominal Closed-Loop Analysis [HINFSYN]\n')
fprintf('=================================================\n');

% Unweight the closedloop plant
CLnom = blkdiag(1/Wu,1/WE)*CLnom*blkdiag(1/Wd,1/Wn);

% Compute L2toE norm for th1 and th2 states from d1 and d2
[gtvn,dWc] = tvnorm(CLnom(3:4,1:2),2,tvnopt);
fprintf(' Euclidean norm bound for theta1 and theta2 due to disturbance d1 and d2:');
disp(gtvn(2));

% Normalized Disturbances d1 and d2
dCL = dScl*dWc/tvnorm(dWc);

% Linear Simulation with Worst-Case Disturbance
[yLin,xLin] = tvlsim(CLnom(3:4,1:2),dCL,tvopt);
yf = tvsubs(yLin,yLin.Time(end));
gLin = norm(yf) / tvnorm(dCL);
fprintf(' Nominal CL Gain Linear Sim = %4.4f\n',gLin);
etaLin = etabar;
etaLin.InterpolationMethod = 'Linear';
etaLin = etaLin +  evalt(xLin(1:4),Tgrid);

% Nonlinear Simulation
d = dCL;
KOpFb = Knom;
NoisePowerSim = [0;0];
model = 'TwoLinkRobotCL_OpFb';
sim(model,[T0 Tf]);
yfNL = x(end,1:2);
gNL = norm(yfNL) / tvnorm(dCL);
fprintf(' Nominal CL Gain NonLinear Sim = %4.4f\n',gNL);
etaNL = tvmat(eta',tsim);
etaNL = evalt(etaNL,Tgrid);

figure;clf;
plot(etabar(1),etabar(2),'k-.', etaNL(1),etaNL(2),'r', ...
    etaLin(1),etaLin(2),'b--');
xlabel('\theta_1 (rad)');
ylabel('\theta_2 (rad)');
legend('Trim','NL','Lin');
title(sprintf('Simulations with ||dCL|| = %.2f',dScl));

figure;clf;
tvplot(dCL);
grid on; box on;
legend('d_1','d_2');
xlabel('Time (s)');
ylabel('Disturbance Magnitude');
title(sprintf('Worst-Case Disturbance ||dCL|| = %.2f',dScl));

%% Save Data
save(mfilename,'Knom','Gnom','gnom','nominfo');