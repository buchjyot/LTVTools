%% Load Example Data
% Second order LTI system
G = ss(tf(1,[1,0.6,1]));

% Time-varying representation on long enough horizon
T0 = 0;
Tf = 300;
t  = T0:0.01:Tf;
Gt = tvss(G,[T0,Tf]);

%% Options
OdeSolver = 'ode45';
tOpt   = tvodeOptions('OdeSolver',OdeSolver);
pOpt   = poweritOptions('Display','on','StopTol',1e-2,'OdeSolver',OdeSolver,...
        'StoreAllIter',true);
tvnOpt = tvnormOptions('OdeSolver',OdeSolver);

%% Infinite Horizon Results
% Hinfnorm and corresponsing disturbance
[g1,w] = hinfnorm(G);
figure; bodemag(G);grid on;

% Plot disturbance
figure;
d = tvmat(sin(w*t),t);
tvplot(d,'LineWidth',2);grid on;box on;
ylabel('Disturbance Input');ylim([-2,2])

% Verify L2 Gain
dn = d/tvnorm(d);
Y = tvlsim(Gt,dn,tOpt);
g1sim = tvnorm(Y); % must be very close to g1

% Display
fprintf('### Infinite Horizon:\n');
fprintf(' HinfNorm: %.3f, Finite Horizon Sim: %.3f\n\n',g1,g1sim);

%% Power Itertaions
fprintf('### Power Iterations:\n');
[gLB1,dwc1,info1] = powerit(Gt,[T0,Tf],pOpt);
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);
analyzeResults(info1);

% Plot disturbance
figure;
tvplot(dwc1,'LineWidth',2);grid on;box on;
ylabel('Disturbance Input');

%% Riccati Approach
fprintf('### Riccati Approach:\n');
[gbnd,dwc2,info2] = tvnorm(Gt,tvnOpt);
fprintf(' Lower Bound: %.3f, Upper Bound: %.3f, Computational Time: %.3f seconds\n\n',gbnd(1),gbnd(2),info2.TotalTime);

% Plot disturbance
figure;
tvplot(dwc2,'LineWidth',2);grid on;box on;
ylabel('Disturbance Input')