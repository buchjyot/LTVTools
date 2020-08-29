%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Horizon
T0 = 0;
Tf = 5;

% SISO Plant
Nx = 2;
Ny = 1;
Nu = 1;
G = ss(tf(1,[1 0.2 1]));

% Time-varying representation
Gt = tvss(G,[T0,Tf]);

% Options
NE = 0;
pOpt = poweritOptions('Display','on','OdeSolver','ode45','StoreAllIter',true);
pSpec = poweritSignalSpec('NE',NE);

%% MATLAB Power Itertaions
% MATLAB Perform iterations
fprintf('### MATLAB Power Iterations:\n');
[gLB1,dwc1,info1] = powerit(Gt,[T0,Tf],pSpec,pOpt);

% Analyze Results (This will plot all the results only if StoreAllIter is set to true)
fh = analyzeResults(info1);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

%% Riccati Bisections
% LTVTools
fprintf('### Riccati Approach:\n');
tOpt = tvnormOptions('Display','on');
[tvn2,dwc2,info2] = tvnorm(Gt,NE,tOpt);

% Display
fprintf(' Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n\n',tvn2(1),info2.TotalBisections,info2.TotalTime);

%% Plot
figure;clf;hold on;box on;
fh1 = tvplot(dwc1,'b'  ,'LineWidth',2.5);
fh2 = tvplot(dwc2,'r--','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Worst-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1)],{'Power Iteration','Riccati Bisection'});