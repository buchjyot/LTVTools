%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Example Data
LTISysEx

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

%% MATLAB Power Iterations
% MATLAB Perform iterations
fprintf('### MATLAB Power Iterations:\n');
pOpt = tvpoweritOptions('Display','on','StepSize','Auto');
[gLB1,dwc1,info1] = tvpowerit(Gt,3,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

%% Simulink Power Itertaions
fprintf('### Simulink Power Iterations:\n');
sOpt = slpoweritOptions('Display','on','StepSize','Auto','Objective','L2toE');
[gLB2,dwc2,info2] = slpowerit('LTISys',sOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB2,info2.TotalIter,info2.TotalTime);

%% Plot
figure(1);clf;hold on;box on;
fh1 = tvplot(dwc1,'b'  ,'LineWidth',2.5);
fh2 = tvplot(-dwc2,'r--','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1)],{'MATLAB Power Iteration','Simulink Power Iteration'});