%% TVGRAM SISO Example
% Use tvgram command to compute time-varying controllability and
% observability gramian for SISO LTI 2 by 2 system

% Specify Time
T0 = 0;
Tf = 5;
t = linspace(T0,Tf,1000)';

% Specify tvodeOptions
tvgramopt = tvodeOptions('OdeSolver','ode45',...
    'OdeOptions',odeset('RelTol',1e-4,'AbsTol',1e-6));

%% Stable System

% State Matrices
A = [-4 1; 1 -7];
B = [1;1];
C = eye(2);
D = [0;0];

% State-Space Objects
Gss = ss(A,B,C,D);
G = tvss(Gss);
G = evalt(G,t);

% GRAM
Wc_stableih = gram(Gss,'c')
Wo_stableih = gram(Gss,'O')

% TVGRAM
Wc_stablefh = tvgram(G,'c');
Wo_stablefh = tvgram(G,'O');
Wc_stablefh_end = tvsubs(Wc_stablefh,Tf);
Wo_stablefh_start = tvsubs(Wo_stablefh,T0);

% Plotting
figure;
tvplot(Wc_stablefh,'LineWidth',2);
title('Controllability Gramian (Stable System)');
legend('Wc_{11}','Wc_{12}','Wc_{21}','Wc_{22}');

figure;
tvplot(Wo_stablefh,'LineWidth',2);
title('Observability Gramian (Stable System)');
legend('Wo_{11}','Wo_{12}','Wo_{21}','Wo_{22}');