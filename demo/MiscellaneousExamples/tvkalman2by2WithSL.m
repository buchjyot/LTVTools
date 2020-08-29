%% TVKALMAN SISO Example
% tvkalman with 2by2 plant compare MATLAB and Simulink Results

%% System Specifications
T0 = 0;
Tf = 15;
t = T0:0.01:Tf;
Tspan = [T0,Tf];

% TVSS
A1 = tvmat(-1-2*t,t);
A2 = tvmat(-t,t);
% A = diag([A1 A2]);
A = [A1 -1;-2 A2];
B = [1;1];
C = eye(2);
D = [0;0];
SYS = tvss(A,B,C,D);

% Noise Influence
% Process Noise Influence
Bw = [1;1];
Bv = [0;0];

% Measurement Noise Influence
Dw = [0;0];
Dv = [1;1];

% Process Noise Covariance
QnDesign = diag([1,1]);

% Measurement Noise Covariance
RnDesign = diag([1,1]);

% Initial condition covariance
P0 = diag([1,1]);

%% Time-Varying Kalman Filter
SYSKF = tvss(A,[B Bw],C,[D Dw]);
[Kfh,Lfh,Pfh] = tvkalman(SYSKF,QnDesign,RnDesign,P0,Tspan);

%% Evaluate Kalman Filter Performance
Aaug = [...
    A      zeros(size(A-Lfh*C));...
    Lfh*C    A-Lfh*C;...
    ];
w_zeros = B*0;
w_ones = [1;1];
v_zeros = D*0;
v_ones = [1;1];
Baug = [...
    B w_ones  v_zeros;...
    B w_zeros v_zeros];
Caug = [C C];
Daug = [D w_zeros v_ones];
aug_sys = tvss(Aaug,Baug,Caug,Daug);

%% Simulation Cases
% Case 1: omega = 1; Qn = 1; Rn = 1;
% Case 2: omega = 10; Qn = 0.1; Rn = 0.1;
% Case 3: omega = 1; Qn = 0; Rn = 0;

% Load Simulink Model
model = 'KalmanFilterTest';
load_system(model);

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
    x0 = [0;0];
    xhat0 = [0.5;-0.5];
    X0 = [x0;xhat0];
    
    % Simulate System
    [~,X] = tvlsim(aug_sys,U,U.Time,X0);
    X = evalt(X,t);
    
    % Plot Results
    figure;h=plot(t,X(1:2),'.',t,X(3:4),'-');
    title('Kalman Filter State Estimate Results (ML)');xlabel('Time(s)');
    legend('$x_1$','$x_2$','$\hat{x}_1$','$\hat{x}_2$','Interpreter','Latex');
    % print(sprintf('twoBytwoLTVKF_ML_FIG_%d',i),'-dpdf','-bestfit');
    
    % Simulate Model
    sim(model,t);
    figure;
    h2 = plot(tSim,XSim(:,3:4),'.',tSim,XSim(:,1:2),'-');
    title('Kalman Filter State Estimate Results (SL)');xlabel('Time(s)');
    legend('$x_1$','$x_2$','$\hat{x}_1$','$\hat{x}_2$','Interpreter','Latex');
    % print(sprintf('twoBytwoLTVKF_SL_FIG_%d',i),'-dpdf','-bestfit');
end