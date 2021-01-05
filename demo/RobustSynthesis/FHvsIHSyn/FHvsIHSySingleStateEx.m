% NOTE: Shipping MATLAB does not allow us to have AutoScalingOrder =[0, 0]
% for dksynOptions, but internal code handles this case. Thus, to use this
% example, user need to first start MATLAB as run in administator mode,
% make appropriate minor change to the file that constructs dksynOptions
% and then run this example.

%% Measurement-Feedback Uncertain Interconnection
rng(0);
P  = ss(0.5,1,1,0);

% Independent IO Dimentions
Nx = order(P);
[Ny,Nu] = size(P);

% Dependent IO dimentions
% By default penalize both the states as output
Nex = Nx;
Neu = Nu;
Ne  = Nex + Neu;
Nv = 1;
Nw = 1;

% Choose outputs that are to be penalized in terminal Euclidean sense
NE = Nex;

% Input Design Weights
WdScl = 0.1;
WnScl = 0.01;

% Output Design Weights
WeuScl = 0.2;
WexScl = 1;

% Uncertainty Design Weights
DelInScl = 1;
DelOutScl = 1;

% Promote to state-space objects
Wd = ss(WdScl*eye(Nu));
Wn = ss(WnScl*eye(Ny));
Weu = ss(WeuScl*eye(Nu));
Wex = ss(WexScl*eye(Nex));
WDelIn = ss(DelInScl*eye(Nv));
WDelOut = ss(DelOutScl*eye(Nw));

% System Interconnection
systemnames = 'P Wd Wn Weu Wex WDelIn WDelOut'; %#ok<*NASGU>
inputvar = '[w; d; n; u]';
outputvar = '[WDelIn; Weu; Wex; P(1)+Wn]';
input_to_P = '[WDelOut+u+Wd]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
input_to_Weu = '[u]';
input_to_Wex = '[P]';
input_to_WDelIn = '[Wd+u]';
input_to_WDelOut = '[w]';
cleanupsysic = 'yes';
Gunc = sysic;

% Performance
NE = 0;
[NY,NU] = size(Gunc);

% Uncertainty
beta = 0.6;
Nw = 1;
Nv = 1;
Del = ultidyn('Del',[Nv Nw],'Bound',beta);
Delta = udyn('Delta',[Nw Nv],'UserData',[0,-10,1]);

% Options
opt = dksynOptions;
opt.AutoScalingOrder = [0,0];
opt.DisplayWhileAutoIter = 'off';

% Options for finite horizon synthesis
Display     = 'off';
OdeSolver   = 'ode23s';
Bounds      = [0 5];
RelTol      = 5e-3;
AbsTol      = 1e-4;

tvopt   = tvodeOptions('OdeSolver',OdeSolver);
tvhopt  = tvhinfsynOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvnopt  = tvnormOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvwcopt = tvwcOptions('RDEOptions',tvnopt,'Display','on','MaxIter',15,'Nlmi',20,'Nsp',20,'StopTol',5e-3);
tvropt = tvrobsynOptions('MaxIter',10,'SynthesisOptions',tvhopt,'AnalysisOptions',tvwcopt,...
    'Display','on','DebugMode',false,'StopWhenWithinTol',true);

%% Save
save(mfilename);