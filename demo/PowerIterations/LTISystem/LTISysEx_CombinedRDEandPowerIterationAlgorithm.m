%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Horizon
T0 = 0;
Tf = 60;

% SISO Plant
Nx = 2;
Ny = 1;
Nu = 1;
G = ss(tf(1,[1 2*0.4 1]));

% Time-varying representation
Gt = tvss(G,[T0,Tf]);

% Options
NE = 0;

%% Test Combined Algorithm
% LTVTools
tOpt = tvnormOptions('Display','on','OdeSolver','ode45');
[tvn,dwc,info] = tvnorm(Gt,NE,tOpt);

% Display
fprintf(' Bounds: [%.3f,%.3f], Total RDE Count: %d, Computational Time: %.3f seconds\n\n',tvn(1),tvn(2),info.RDEcnt,info.TotalTime);