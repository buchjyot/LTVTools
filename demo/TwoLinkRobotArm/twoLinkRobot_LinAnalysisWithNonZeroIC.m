%% twoLinkRobot_LinAnalysisWithNonZeroIC
% (Interconnected with all 100 Deltas)
%
% This example uses the method described in the following paper:
%
% Schweidel, Buch, Seiler, Arcak, Computing Worst-Case Disturbances for
% Finite-Horizon Linear Time-Varying Approximations of Uncertain Systems,
% Submitted to LCSS/ACC 2021

%% Load Data
load('twoLinkRobot_BuildLTVModel.mat')
load('twoLinkRobot_LQR.mat')
DeltaStruct = load('RandomDelta.mat','DeltaSim');
DeltaArray = DeltaStruct.DeltaSim;

% Evaluate on same time grid
[G,Klqr] = evalt(G,Klqr,union(G.Time,Klqr.Time));

%% System with Affine Term as Input
systemnames = 'G Klqr';
inputvar = '[w; d(2); vbar]';
outputvar = '[vbar+d(2)-Klqr(2); G(1:2)]';
input_to_G = '[d(1)-Klqr(1); d(2)-Klqr(2)+w]';
input_to_Klqr = '[G]';
cleanupsysic = 'yes';
Tnom_wAff = sysic;

% Scale Plant
Lscl = diag([sqrt(DelNorm) 1 1]);
Rscl = diag([sqrt(DelNorm) 1 1 1]);
Tnom_wAff_scaled = Lscl*Tnom_wAff*Rscl;
tVec = Tnom_wAff_scaled.Time;

%% Sample Uncertainty Delta to Approximate Lower Bound on Worst-Case Gain

% Number of Euclidean Outputs
NE = 2;

% Number of Delta
nD = length(DeltaArray);

% Size of disturbance
alphaNorm = 0.5;

% Memory Allocation
ENormBndNLfromNonzero  = zeros(nD,1);
ENormBndNLfromZero     = zeros(nD,1);
eMaxNonzero            = zeros(nD,1);
eMaxZero               = zeros(nD,1);
norm_eT2a              = zeros(nD,1);
norm_eT2b              = zeros(nD,1);
norm_eT3a              = zeros(nD,1);
norm_eT3b              = zeros(nD,1);

% nonlinear, continuous-time
model = 'TwoLinkRobotCL_Sfb';
load_system(model);

% Trim Input
taubar.InterpolationMethod = 'Linear';
vbar = evalt(taubar(2),tVec);

% Options
tvnOpt = tvnormOptions('Display','off');

%% Main For Loop
for i = 1:nD
    
    %% Display Iteration Count
    disp(i);
    
    %% Augmented system
    % Normalize Uncertainty
    Delta = DeltaArray{i} / norm(DeltaArray{i},inf);
    
    % Augmented System
    sysAug = augsyswithunc(Tnom_wAff_scaled,Delta,vbar,NE);
    
    %% Analysis with Zero IC
    [gZ,dZ] = tvnorm(sysAug,NE,tvnOpt);
    
    % Returned output is unit norm disturbance and related bound
    dStarZero   = dZ*alphaNorm;
    eMaxZero(i) = gZ*alphaNorm;
    
    %% Analysis with Non-zero IC
    [eMaxNonzero(i),dStarNonzero] = tvnormaff(sysAug,alphaNorm,NE,tvnOpt);
    
    %% Simulate Linear Model
    
    % % With zero IC
    % [e2a, x2a] = tvlsim(sysAug,dStarZero,[T0,Tf]);
    % [e2b, x2b] = tvlsim(sysAug,dStarNonzero,[T0,Tf]);
    %
    % norm_eT2a(i) = norm(e2a.Data(:,:,end));
    % norm_eT2b(i) = norm(e2b.Data(:,:,end));
    %
    % % With nonzero IC
    % xNonzero = sysAug.UserData;
    % [e3a, x3a] = tvlsim(sysAug,dStarZero,[T0,Tf],xNonzero);
    % [e3b, x3b] = tvlsim(sysAug,dStarNonzero,[T0,Tf],xNonzero);
    %
    % norm_eT3a(i) = norm(e3a.Data(:,:,end));
    % norm_eT3b(i) = norm(e3b.Data(:,:,end));
    
    %% Simulate Nonlinear Model
    
    % % Simulink model settings
    % KSfb = -Klqr;
    % Delta = ss(sqrt(DelNorm)*Delta*sqrt(DelNorm)); % how is specific realization chosen?
    % DeltaSim = Delta;
    %
    % % (1a) dStarZero
    % d = dStarZero;
    % sim(model)
    % ENormBndNLfromZero(i) = norm(x(end,1:2));
    %
    % % (1b) dStarNonzero
    % d = dStarNonzero;
    % sim(model)
    % ENormBndNLfromNonzero(i) = norm(x(end,1:2));
end

%% Save Data
save('twoLinkRobot_LinAnalysisWithNonZeroIC',...
    'nD','eMaxZero','eMaxNonzero');%,'ENormBndNLfromZero','ENormBndNLfromNonzero',...
    %'norm_eT2a','norm_eT2b','norm_eT3a','norm_eT3b');
return;

%% Plot
load('twoLinkRobot_LinAnalysisWithNonZeroIC.mat');

% Plot 1: Linear Sim
[nonzeroSortL,ind] = sort(eMaxNonzero);
zeroSortL = eMaxZero(ind);
figure;clf;hold on; grid on;box on;
plot(1:nD,nonzeroSortL.^2, 'rd',1:nD,zeroSortL.^2, 'bo');
title('Euclidean Norm Bound from Linear Analysis');
ylabel('$R_\alpha$','fontsize',13,'Interpreter','latex');xlabel('\Delta index (sorted)');
legend({'$R_\alpha^i$','$R_\alpha^{0,i}$'},'fontsize',13,'Interpreter','latex','Location','northwest');

% % Plot 2: Nonlinear Sim
% [nonzeroSortNL,ind] = sort(ENormBndNLfromNonzero);
% zeroSortNL = ENormBndNLfromZero(ind);
% figure;clf;hold on; grid on;box on;
% plot(1:nD,nonzeroSortNL, 'rd',1:nD,zeroSortNL, 'bo');
% title('Euclidean Norm Bound from Nonlinear Simulation');
% xlabel('\Delta index (sorted)');
% ylabel('$\|e(T)\|_2$','Interpreter','latex');
% legend('NLSim with d_\alpha','NLSim with d_\alpha^0','Location','northwest');