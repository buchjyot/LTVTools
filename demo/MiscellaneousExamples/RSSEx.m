%% RSSEx
%
% This file contains the example data for the 4 state Random LTI problem to
% be analysed on Finite Horizon

%% Specify Options
Display = 'on';
RelTol = 1e-4;
AbsTol = 1e-4;
Bounds = [0 5];
OdeSolver = 'ode23s';

% Options
tvhopt = tvhinfsynOptions('Bounds',Bounds,'Display',Display,...
    'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);

tvnopt = tvnormOptions('Bounds',Bounds,'Display',Display,'OdeSolver',OdeSolver,...
    'RelTol',RelTol,'AbsTol',AbsTol);

tvopt = tvodeOptions('OdeSolver',OdeSolver);

% Time Horizon
T0 = 0;
Ts = 0.1;
Tf = 1;

%% System Dimentions
Nx = 4;
Ny = 2;
Nu = 1;
Nd1 = 1;
Nex = 4;
Neu = 1;
Ne = Nex+Neu;

%% Random LTI State-Space
% xdot = A*x + Bu*u
% y    = C*x + Du*u

% Choose Random Plant Model
rng('shuffle');
P = rss(Nx,Ny,Nu);
[A,B,C,D] = ssdata(P);
A = round(A*10)/10;
B = round(B*10)/10;
C = round(C*10)/10;
D = 0*round(D*10)/10;
P = ss(A,B,C,D);

% State Matrices
% A = diag([-4 -3 -2 -1]);
% B = ones(4,1);
% C = [zeros(2) eye(2)];
% D = zeros(2,1);
% P = ss(A,B,C,D);

%% Disturbance and Noise Weights
Wd = ss(0.1*eye(Nd1));
Wn = ss(0.01*eye(Ny));

%% Full-Information Problem
% Add disturbances

% SISIC
systemnames = 'P Wd'; %#ok<*NASGU>
inputvar = '[d; u]';
outputvar = '[P]';
input_to_P = '[u+Wd]';
input_to_Wd = '[d]';
cleanupsysic = 'yes';
G = sysic;

% Choose Outputs
NU = size(G,2);
Cfi = [zeros(Nu,Nx);eye(Nx)];
Dfi = [zeros(Nu,NU-Nu) eye(Nu);zeros(Nx,NU)];

% Full-Information Problem
Gfi = ss(G.A,G.B,Cfi,Dfi);

%% Output-Feedback Problem
% Add process noise p(t) and measurement noise n(t)

% SYSIC
systemnames = 'P Wd Wn';
inputvar = '[d; n(2); u]';
outputvar = '[P+Wn]';
input_to_P = '[u+Wd]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
cleanupsysic = 'yes';
G = sysic;

% Choose Outputs
NU = size(G,2);
C1 = [zeros(Nu,Nx);eye(Nx)];
Cnom = [C1;G.C];
Dnom = [zeros(Nu,NU-Nu) eye(Nu);zeros(Nx,NU);G.D];

% Nominal Plant with all the disturbances
Gnom = ss(G.A,G.B,Cnom,Dnom);

%% Uncertain Plant Interconnection
% w is uncertain input,
% d is disturbance input
% n is measurement noise
% u is control input
% v is uncertain output

% SYSIC
systemnames = 'P Wd Wn';
inputvar = '[w; d; n(2); u]';
outputvar = '[u+Wd;P+Wn]';
input_to_P = '[u+Wd+w]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
cleanupsysic = 'yes';
G = sysic;

% IO dimentions
Nw = 1;
Nv = 1;

% Choose Outputs
NU = size(G,2);
Cunc = [G.C(Nv,:);zeros(Nu,Nx);eye(Nx);G.C(end-Ny+1:end,:)];
Dunc = [G.D(Nv,:);zeros(Nu,NU-Nu) eye(Nu);zeros(Nx,NU);G.D(end-Ny+1:end,:)];

% Uncertain Plant with all the disturbances
Gunc = ss(G.A,G.B,Cunc,Dunc);
[NX,NY,NU] = size(Gunc);

%% Save Data to MAT file
save(mfilename);