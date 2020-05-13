%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Example Data
LTISystemEx

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

% Options
pOpt = poweritOptions('Display','on','StepSize','Auto','InitialInput','ones');

%% MATLAB Power Itertaions
% MATLAB Perform iterations
fprintf('### MATLAB Power Iterations:\n');
[gLB1,dwc1,info1] = powerit(Gt,pOpt);

% Display
fprintf('L2 Gain Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

%% Simulink Power Itertaions
% Simulink Perform iterations
fprintf('### SIMULINK Power Iterations:\n');
model = 'LTISys';
[gLB2,dwc2,info2] = powerit(model,pOpt);

% Display
fprintf('L2 Gain Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB2,info2.TotalIter,info2.TotalTime);

%% Riccati Bisections
% LTVTools
[tvn,dwc3,info3] = tvnorm(Gt);

% Display
fprintf('### Riccati Approach:\n');
fprintf('L2 Gain Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n',tvn(1),info3.TotalBisections,info3.TotalTime);

%% Plot
figure(1);clf;hold on;box on;
fh1 = tvplot(dwc1,'b'  ,'LineWidth',2.5);
fh2 = tvplot(dwc2,'r--','LineWidth',2.5);
fh3 = tvplot(-dwc3,'g-.','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1) fh3(1)],{'MATLAB Power Iteration','Simulink Power Iteration','Riccati Bisection'});