%% TVKALMAN SISO Example
% Use tvkalman to design an LTV Kalman Filter for SISO one by one LTV plant

%% System Specification
T0 = 0;
Tf = 15;
t = T0:0.01:Tf;
Tspan = [T0,Tf];
Adata = -5 - 2*t;

% TVSS
A = tvmat(Adata,t);
B = 1;
C = 1;
D = 0;
SYS = tvss(A,B,C,D);

% Noise Influence
% Process Noise Influence
Bw = 1;
Bv = 0;

% Measurement Noise Influence
Dw = 0;
Dv = 1;

% Process Noise Covariance
QnDesign = 1;

% Measurement Noise Covariance
RnDesign = 1;

% Initial condition covariance
P0 = 1;

%% Finite Horizon Kalman Filter
SYSKF = tvss(A,[B Bw],C,[D Dw]);
[Kfh,Lfh,Pih] = tvkalman(SYSKF,QnDesign,RnDesign,P0,Tspan);

%% Evaluate Kalman Filter Performance
Aaug = [...
    A      zeros(size(A-Lfh*C));...
    Lfh*C    A-Lfh*C;...
    ];
Baug = [...
    B 1  0;...
    B 0 0];
Caug = [C C*0];
Daug = [D 0 1];
aug_sys = tvss(Aaug,Baug,Caug,Daug);

% Simulation Cases
% Case 1: omega = 1; Qn = 1; Rn = 1;
% Case 2: omega = 10; Qn = 1; Rn = 1;
% Case 3: omega = 1; Qn = 0; Rn = 0;

for i = 1:3
    switch i
        case 1
            omega = 1; Qn = 1; Rn = 1;
        case 2
            omega = 10; Qn = 0.1; Rn = 0.1;
        case 3
            omega = 1; Qn = 0; Rn = 0;
    end
    
    % Create Input Signals
    n = length(t);
    W = tvmat(sqrt(Qn)*randn(n,1),t);
    V = tvmat(sqrt(Rn)*randn(n,1),t);
    
    % Control input is 0, initial condition response
    u = tvmat(sin(omega*t),t);
    U = [u;W;V];
    
    % Augmented Plant initial condition
    X0 = [0;0.3];
    
    % Simulate System
    [~,X] = tvlsim(aug_sys,U,X0);
    X = evalt(X,t);
    
    % Plot Results
    figure;plot(t,X(1),'.',t,X(2),'-');
    title('Kalman Filter State Estimate Results');
    legend('x','xhat');xlabel('Time(s)');
    % print(sprintf('oneByoneLTVKF_ML_FIG_%d',i),'-dpdf','-bestfit');
end