% NOTE: Shipping MATLAB does not allow us to have AutoScalingOrder =[0, 0]
% for dksynOptions, but internal code handles this case. Thus, to use this
% example, user need to first start MATLAB as run in administator mode,
% make appropriate minor change to the file that constructs dksynOptions
% and then run this example.

%% Measurement-Feedback Uncertain Interconnection
rng(0);

% Mass-Spring-Damper System
% m := mass
% b := damping coefficient
% k := spring constant
msdEx1 = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
P  = msdEx1(1,0.8,1);

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

% Transform plant to have all the states (Nx) as outputs
P = ss(P.A,P.B,eye(Nx),0);

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

%% Perform Infinite Horizon Robust Synthesis
Gu = lft(Del,Gunc);
Gscl = Gu;

gUpp = 2;
gLow = 0.01;
gTry = (gLow+gUpp)/2;
Tol = 1e-4;

while true
    % Scale performance channel until Robust Performance from DKsyn is 1
    Gscl = Gu*blkdiag(eye(NU-Nw-Ny)/gTry,eye(Ny));
    
    % Run DK Iterations
    [K,CL,GAM,DKINFO] = dksyn(Gscl,Ny,Nu,opt);
    gstr = dksynperf(CL);
    gRP  = gstr.UpperBound;
    
    % Bisect until RobustPerformance ~= 1
    fprintf('gUpp = %.4f, gLow = %.4f, gTry = %.4f, gRP = %.4f\n',gUpp,gLow,gTry,gRP);
    if abs(gRP-1) <= Tol
        break;
    end
    if gRP > 1
        gLow = gTry;
    else
        gUpp = gTry;
    end
    gTry = (gLow+gUpp)/2;
end
gWCGain = gTry;
fprintf('LTI Synthesis worst-case gain : %.4f\n',gWCGain);

%% Finite Horizon L2toL2 Synthesis
Tall = 200; %50, 20, 10, 5, 3, 1];
NT = length(Tall);
Krob = cell(NT,1);
CLrob = cell(NT,1);
grob = zeros(NT,1);

Display     = 'off';
OdeSolver   = 'ode23s';
Bounds      = [0 5];
RelTol      = 5e-3;
AbsTol      = 1e-4;

tvopt   = tvodeOptions('OdeSolver',OdeSolver);
tvhopt  = tvhinfsynOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvnopt  = tvnormOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvwcopt = tvwcOptions('RDEOptions',tvnopt,'Display','on','MaxIter',15,'Nlmi',20,'Nsp',20,'StopTol',5e-3);
tvropt = tvrobsynOptions('MaxIter',20,'SynthesisOptions',tvhopt,'AnalysisOptions',tvwcopt,...
    'Display','on','DebugMode',false,'StopWhenWithinTol',true);

GuncSCL = blkdiag(sqrt(beta),eye(NY-Nv))*Gunc*blkdiag(sqrt(beta),eye(NU-Nw));
for i = 1:NT
    Gtv = tvss(GuncSCL,linspace(0,Tall(i),Tall(i)*2));
    [Krob{i},CLrob{i},grob(i)] = tvrobsyn(Gtv,Delta,Ny,Nu,NE,tvropt);
end

%% Save
save(mfilename);