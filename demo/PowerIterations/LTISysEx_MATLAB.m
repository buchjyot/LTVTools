%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Example Data
LTISysEx

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

% Number of Euclidean Outputs
NE = 0;

%% MATLAB Power Itertaions
% MATLAB Perform iterations
fprintf('### MATLAB Power Iterations:\n');
pOpt = tvpoweritOptions('Display','on','StepSize','Auto');
[gLB1,dwc1,info1] = tvpowerit(Gt,NE,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

%% Riccati Bisections
% LTVTools
fprintf('### Riccati Approach:\n');
tOpt = tvnormOptions('Display','on');
[tvn2,dwc2,info2] = tvnorm(Gt,NE,tOpt);

% Display
fprintf(' Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n\n',tvn2(1),info2.TotalBisections,info2.TotalTime);

%% Best of Both
% LTVTools
fprintf('### Best of Both:\n');
tOpt = tvnormOptions('Display','on');
[tvn3,dwc3,info3] = tvnorm1(Gt,NE,tOpt);

% Display
fprintf(' Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n\n',tvn3(1),info3.TotalBisections,info3.TotalTime);

%% Plot
figure(1);clf;hold on;box on;
fh1 = tvplot(dwc1,'b'  ,'LineWidth',2.5);
fh2 = tvplot(dwc2,'r--','LineWidth',2.5);
fh3 = tvplot(dwc3,'g-.','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1) fh3(1)],{'Power Iteration','Riccati Bisection','Combined'});