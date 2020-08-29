%% TwoLinkRobot_PowerIt
% This example uses the two-link robot arm example for the power iterations
% to consturct worst-case disturbance using the nonlinear and LTV model. It
% compares the obtained disturbance with LTV analysis using Riccati
% approach. We use the MAT file data logged in the TwoLinkRobot demo folder

%% Load Data
% These files are in the TwoLinkRobot demo folder.
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_LQR.mat');
KSfb = -Klqr;
NE = 2;

%% LTV Analysis
% Nominal Closed-Loop Loop L2 to E Gain
% Expeced 0.05509 (MS Thesis)
[gCLNom,dCL,infoCL] = tvnorm(Tlqr(1:2,:),NE);
fprintf(' Bounds on Nominal CL Gain = [%4.4f, %4.4f] \n',...
    gCLNom(1),gCLNom(2));
dNorm = 1;
dCL = dNorm*dCL/tvnorm(dCL);

%% Power Iterations
model = 'TwoLinkRobotNLEx';
pOpt  = poweritOptions('Display','on','StopTol',1e-2,'Display','on','StoreAllIter',true);
pSpec = poweritSignalSpec('NE',NE);

disturbanceInitCase = 2;
tgrid = dCL.Time;
Nt = length(tgrid);
switch disturbanceInitCase
    
    case 1
        % It takes many iterations (specifically 202 iterations in LTV case) 
        % for the disturbance to converge to the worst-case disturbance.
        % where as the performance converges within very few (3 or 4) iterations. 
        % Thus, we may initialize the power iteration algorithm with the
        % disturbance that we got from LTV analysis to speed up the
        % iterations.
        U1 = dCL;
        
    case 2
        % Let the algorithm take as many as iterations it take. 
        U1 = tvmat(randn(2,1,Nt),tgrid);
        
    case 3
        % Specify defaults
        U1 = [];
        pSpec.InitialInput = 'ones';
end
        
% LTV vs Nonlinear power iterations
poweritCase = 1;
switch poweritCase
    case 1
        % LTV 
        [T0,Tf] = getHorizon(U1);
        [glb,dwc,info] = powerit(Tlqr(1:2,:),[T0,Tf],U1,pSpec,pOpt);

    case 2
        % Nonlinear 
        
        % NOTE: Simulink linearization perturbs the nonlinear ode solutions
        % along the trajectory to compute (A,B,C,D) matrices. We use the
        % TVMAT block inside the simulink model. We throw warnings if we
        % try to evaluate outside of the regime. Thus, for now the horizon
        % needs to be something less than the boundry conditions in order
        % to allow room for the numerical perturbation.
        U1 = evalt(U1,[0.01 0.01:0.01:4.99 4.99]);
        [T0,Tf] = getHorizon(U1);
        [glb,dwc,info] = powerit(model,[T0,Tf],U1,pSpec,pOpt);
        
end
fprintf(' Power Iteration Lower Bound = [%4.4f] \n',glb);

%% Plot
figure(1);clf;hold on;box on;grid on;
fh1 = tvplot(dwc(1),'b',dwc(2),'r','LineWidth',2);
fh2 = tvplot(dCL(1),'b--',dCL(2),'r--','LineWidth',2);
xlabel('Time (sec)');
ylabel('Worst-Case Disturbance Signal');
title('Solid: Power Iterations, Dashed: LTV Analysis');