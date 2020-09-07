%% TVGRAM SISO Example
% Use tvgram command to compute time-varying controllability and
% observability gramian for SISO system

% Specify Time
T0 = 0;
Tf = 5;
t = linspace(T0,Tf,100)';

% Specify tvodeOptions
tvgramopt = tvodeOptions('OdeSolver','ode45',...
    'OdeOptions',odeset('RelTol',1e-4,'AbsTol',1e-6));

%% Stable System

% State Matrices
AData = -5 - 2*t; % Negative eigen values
A = tvmat(AData,t);
B = 2;
C = 1;
D = 0;

% State-Space Objects
G = tvss(A,B,C,D);

% TVGRAM
Wc_stable = tvgram(G,'c',tvgramopt);
Wo_stable = tvgram(G,'O',tvgramopt);

% Plotting
figure;
plot(Wc_stable.Time,Wc_stable,'LineWidth',2);
title('Controllability Gramian (Stable System)');

figure;
plot(Wo_stable.Time,Wo_stable,'LineWidth',2);
title('Observability Gramian (Stable System)');

%% Unstable System

% State Matrices
AData = 5 + 2*t; % Positive values
A = tvmat(AData,t);

% TVSS
G = tvss(A,B,C,D);

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