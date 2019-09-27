%% TVGRAM SISO Example
% Use tvgram command to compute time-varying controllability and
% observability gramian for SISO LTI system

% Specify Time
T0 = 0;
Tf = 5;
t = linspace(T0,Tf,1000)';

% Specify tvodeOptions
tvgramopt = tvodeOptions('OdeSolver','ode45',...
    'OdeOptions',odeset('RelTol',1e-4,'AbsTol',1e-6));

%% Stable System

% State Matrices
A = -4;
B = 2;
C = 1;
D = 0;

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
plot(Wc_stablefh.Time,Wc_stablefh,'LineWidth',2);
title('Controllability Gramian (Stable System)');

figure;
plot(Wo_stablefh.Time,Wo_stablefh,'LineWidth',2);
title('Observability Gramian (Stable System)');

%% Unstable System

% State Matrices
A = 5; % Positive values

% State-Space Objects
Gss = ss(A,B,C,D);
G = tvss(Gss);
G = evalt(G,t);

% GRAM
% Wc_unstableih = gram(Gss,'c'); % Error
% Wo_unstableih = gram(Gss,'O');

% TVGRAM
Wc_unstable = tvgram(G,'c',tvgramopt);
Wo_unstable = tvgram(G,'O',tvgramopt);

% Plotting
figure;
plot(Wc_unstable.Time,Wc_unstable,'LineWidth',2);
title('Controllability Gramian (Unstable System)');

figure;
plot(Wo_unstable.Time,Wo_unstable,'LineWidth',2);
title('Observability Gramian (Unstable System)');