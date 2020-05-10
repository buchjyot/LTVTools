%% TwoStateEx
% 2-state numerical (LTI) example on finite-horizon for robust synthesis

%% Problem Data
% Allows user to try different plants
plant_case = 1;

switch plant_case
    case 1
        %% Mass-Spring-Damper System
        % m := mass
        % b := damping coefficient
        % k := spring constant
        msdEx1 = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
        P  = msdEx1(1,0.8,1);
end

% Independent IO Dimentions
Nx = order(P);
[Ny,Nu] = size(P);

% Dependent IO dimentions
% By default penalize both the states as output
Nex = Nx;
Neu = Nu;
Ne  = Nex + Neu;

% Choose outputs that are to be penalized in terminal Euclidean sense
NE = Nex;

% Transform plant to have all the states (Nx) as outputs
P = ss(P.A,P.B,eye(Nx),0);

% Horizon
T0 = 0;
Ts = 0.1;
Tf = 3;

% Uncertain IO dimentions
Nv = 1;
Nw = 1;

%% Measurement-Feedback Uncertain Interconnection
% Input Design Weights
WdScl = 0.1;
WnScl = 0.01;

% Output Design Weights
WeuScl = 1;
WexScl = 2;

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
Gs = sysic;

% Finite Horizon Representations
Gunc = evalt(tvss(Gs),T0:Ts:Tf);
Gnom = lft(0,Gunc);

%% State-Feedback Uncertain Interconnection
% System Interconnection
systemnames = 'P Wd Weu Wex WDelIn WDelOut';
inputvar = '[w; d; u]';
outputvar = '[WDelIn; Weu; Wex]';
input_to_P = '[WDelOut+u+Wd]';
input_to_Wd = '[d]';
input_to_Weu = '[u]';
input_to_Wex = '[P]';
input_to_WDelIn = '[Wd+u]';
input_to_WDelOut = '[w]';
cleanupsysic = 'yes';
Gs = sysic;

% Finite Horizon Representation
GuSfb = evalt(tvss(Gs),T0:Ts:Tf);
GnSfb = lft(0,GuSfb);

%% Full-Information Nominal Interconnection
% Full-Information Problem assumes that we have an access to states and
% disturbances, thus add disturbances to the plant and write in terms of
%
% xdot = A(t)x(t) + Bd(t)*d(t) + Bu(t)*u(t)
% e(t) = C(t)x(t) + Dd(t)*d(t) + Du(t)*u(t)

% SISIC
systemnames = 'P Wd Weu Wex'; %#ok<*NASGU>
inputvar = '[d; u]';
outputvar = '[Weu; Wex]';
input_to_P = '[u+Wd]';
input_to_Weu = '[u]';
input_to_Wex = '[P]';
input_to_Wd = '[d]';
cleanupsysic = 'yes';
Gs = sysic;

% Finite Horizon Representation
Gfi = evalt(tvss(Gs),T0:Ts:Tf);

%% Describe Uncertainty using IQCs
% NormBound on Uncertainty Delta
% || Delta || <= beta
Delta = udyn('Delta',[Nw Nv],'UserData',[0,0,0]);
beta = 0.6;

%% Available Options
Display     = 'off';
OdeSolver   = 'ode23s';
Bounds      = [0 5];
RelTol      = 3e-3;
AbsTol      = 1e-4;

tvopt   = tvodeOptions('OdeSolver',OdeSolver);
tvspt   = tvlsimOptions('OdeSolver',OdeSolver);
tvhopt  = tvhinfsynOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvnopt  = tvnormOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvwcopt = tvwcOptions('RDEOptions',tvnopt,'Display','on','MaxIter',10,'Nlmi',25,'StopTol',5e-3);
tvropt = tvrobsynOptions('MaxIter',20,'SynthesisOptions',tvhopt,'AnalysisOptions',tvwcopt,...
    'Display','on','DebugMode',false);

%% Save Data
save(mfilename);