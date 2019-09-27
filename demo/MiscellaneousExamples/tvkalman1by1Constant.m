%% TVKALMAN SISO Example
% Use tvkalman to design an LTV Kalman Filter for SISO one by one constant
% plant

%% System Specification

% Specify Time Grid
T0 = 0;
Tf = 30;
Time = linspace(T0,Tf,1000)';

% State Matrices
a = -3;
A = tvmat(a);
A = evalt(A,Time);
B = 1;
C = 1;
D = 0;
G = 1;
H = 0;

% TVSS
sys_tvss = tvss(A,B,C,D);

% Initial Conditions
P0 = 1;
x0 = -1;

% Process Noise Covariance
Qn = 1;

% Measurement Noise Covariance
Rn = 1;

%% Infinite Horizon Kalman Filter
% For Inf. Horizon kalman filter
%      .
%      x = Ax + Bu + Gw            {State equation}
%      y = Cx + Du + Hw + v        {Measurements}
sys_ss = ss(a,[B G],C,[D H]);

% Use Infinite Horizon kalman filter to compute gains and
% covariances (Ricatti Diffrential Equation solutions)
[Kih,Lih,Pih] = kalman(sys_ss,Qn,Rn);

%% Finite Horizon Kalman Filter
% System Under Test Parameters and time varying system
%      .
%      x = A(t)x(t) + B(t)u(t) + G(t)*w        {State equation}
%      y = C(t)x(t) + D(t)u(t) + H(t)*w + v    {Measurements}
sys_tvss_kf = tvss(A,[B G],C,[D H]);

% Make sure LTV system is Observable
Wo = tvgram(sys_tvss_kf,'o');
figure;
plot(Wo.Time,Wo,'LineWidth',2);
title('Observability Gramian (Stable System)');

% Use TVKALMAN to integrate CDRE
[Kfh,Lfh,Pfh] = tvkalman(sys_tvss_kf,Qn,Rn,P0,[0 Tf]);
[Ae,Be,Ce,De] = ssdata(Kfh);

% Plot the data
figure;
plot(Pfh.Time,Pfh,'LineWidth',2);
title('Solution of Riccati Equation (Error Covariance P)');

figure;
plot(Lfh.Time,Lfh,'LineWidth',2);
title('Kalman Filter Gain Values');

% Evaluate Kalman Filter Performance
[A,B,C,D] = ssdata(sys_tvss);
Aaug = [...
    A      zeros(size(A-Lfh*C));...
    Lfh*C    A-Lfh*C];
Baug = [...
    B 1 0;...
    B 0 0];
Caug = [C C*0];
Daug = [D 0 0];
aug_sys = tvss(Aaug,Baug,Caug,Daug);

% Create Input Signals
t = sys_tvss_kf.Time;
n = length(sys_tvss_kf.Time);
W = tvmat(sqrt(Qn/10)*randn(n,1),t);
V = tvmat(sqrt(Rn/10)*randn(n,1),t);

% Control input is 0, i.e. initial condition response
u = tvmat(zeros(n,1),t);
U = [u;W;V];

% Augmented Plant initial condition
X0 = [x0;x0+4];

% Simulate System
[~,X] = tvlsim(aug_sys,U,X0);
X = evalt(X,t);

% Plot Results
figure;
plot(t,X(1),t,X(2),'LineWidth',2);
title('Kalman Filter State Estimate Results');
legend('x','xhat');