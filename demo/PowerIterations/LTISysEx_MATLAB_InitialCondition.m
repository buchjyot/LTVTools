%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Example Data
LTISysEx

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

% Number of Euclidean Outputs
NE = 0;

% Options
pOpt = tvpoweritOptions('Display','on','StepSize','Auto');

%% Power Itertaions without Initial Condition
fprintf('### MATLAB Power Iterations without Initial Condition:\n');
[gLB1,dwc1,info1] = tvpowerit(Gt,NE,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

%% Power Itertaions with Initial Condition
fprintf('### MATLAB Power Iterations with Initial Condition:\n');
x0 = 5*randn(Nx,1);
[gLB2,dwc2,info2] = tvpowerit(Gt,NE,x0,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB2,info2.TotalIter,info2.TotalTime);

%% Plot
figure(1);clf;hold on;box on;
fh1 = tvplot(dwc1,'b'  ,'LineWidth',2.5);
fh2 = tvplot(dwc2,'r--','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1)],{'Without I.C.','With I.C.'});