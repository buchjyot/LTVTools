%% Kirk Example
% This example uses second order example taken from the Kirk's book to
% design finite horizon LQR

% Reference:
% Kirk, D. E. (2012). Optimal control theory: an introduction. Courier
% Corporation.(Ex 5.2-2) Second Order LTI system

%% System Specifications
% Time Horizon
T0 = 0;
Tf = 15;

% State Matrices
A = [0 1;2 -1];
B = [0; 1];
C = eye(2);
D = 0;

% Noise Influence
Q = [2 0;0 1];
R = 1/2;
F = zeros(2);

% State-Space Objects
sys_ss = ss(A,B,C,D);
sys = tvss(sys_ss);

%% Infinite Horizon LQR
[Kih,Sih] = lqr(sys_ss,Q,R);

%% Finite Horizon LQR
[Kfh,Pfh] = tvlqr(sys,Q,R,F,[T0 Tf]);

%% Plot the data
figure(1);clf;
tvplot(Pfh,'LineWidth',2);grid on;
legend('P11','P12','P21','P22','Location','best');
title('Solution of Riccati Equation');

% NOTE: The P plot should match with Page no. 218, Figure 5-9 (a)
figure(2);clf;
tvplot(Kfh,'LineWidth',2);grid on;
legend('K1','K12','Location','best');
title('State Feedback Gain Values');

% Do simulation and plot state trajectories
x0 = [-4 4];
Acl = A - B*Kfh;
Bcl = sys.B*0;
R = tvmat(zeros(length(Acl.Time),1),Acl.Time); % Track 0
sys_cl = tvss(Acl,Bcl,sys.C,sys.D);
X = tvlsim(sys_cl,R,x0);
Xdata = reshapedata(X);

% Control inputs
Xeval = evalt(X,Kfh.Time);
U = -Kfh*Xeval;
Udata = reshapedata(U);

% Plot Page no. 218, Figure 5-9 (b)
figure(3);clf;
tvplot(X,'--r',U,'b','LineWidth',2); grid on;
legend('X1*','X2*','U*');
title('Optimal Control Input and State Trajectories');