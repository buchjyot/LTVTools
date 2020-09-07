% NOTE: Shipping MATLAB does not allow us to have AutoScalingOrder =[0, 0]
% for dksynOptions, but internal code handles this case. Thus, to use this
% example, user need to first start MATLAB as run in administator mode,
% make appropriate minor change to the file that constructs dksynOptions
% and then run this example.

% Plant
rng(0);
P = ss(0.2,1,1,0);
Ny = 1;
Nu = 1;
Nex = 1;

% Uncertainty
beta = 0.5;
Nw = 1;
Nv = 1;
Del = ultidyn('Del',[Nv Nw],'Bound',beta);
Delta = udyn('Delta',[Nw Nv],'UserData',[0,-10,1]);

% Input Design Weights
WdScl = 0.3;
WnScl = 0.1;

% Output Design Weights
WeuScl = 1;
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
outputvar = '[WDelIn; Weu; Wex; P+Wn]';
input_to_P = '[WDelOut+u+Wd]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
input_to_Weu = '[u]';
input_to_Wex = '[P]';
input_to_WDelIn = '[Wd+u]';
input_to_WDelOut = '[w]';
cleanupsysic = 'yes';
Gunc = sysic;
[NY,NU] = size(Gunc);

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
Tall = [0.5, 0.7, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 7, 10, 20, 50, 100];
NT = length(Tall);
Krob = cell(NT,1);
CLrob = cell(NT,1);
grob = zeros(NT,1);

Display     = 'off';
OdeSolver   = 'ode23s';
Bounds      = [0 5];
RelTol      = 5e-3;
AbsTol      = 5e-4;

tvopt   = tvodeOptions('OdeSolver',OdeSolver);
tvhopt  = tvhinfsynOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvnopt  = tvnormOptions('Bounds',Bounds,'OdeSolver',OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvwcopt = tvwcOptions('RDEOptions',tvnopt,'Display','on','MaxIter',10,'Nlmi',20,'StopTol',5e-3);
tvropt = tvrobsynOptions('MaxIter',15,'SynthesisOptions',tvhopt,'AnalysisOptions',tvwcopt,...
    'Display','on','DebugMode',false,'StopWhenWithinTol',true);

GuncSCL = blkdiag(sqrt(beta),eye(NY-Nv))*Gunc*blkdiag(sqrt(beta),eye(NU-Nw));
for i = 1:NT
    Gtv = tvss(GuncSCL,[0 Tall(i)]);
    [Krob{i},CLrob{i},grob(i)] = tvrobsyn(Gtv,Delta,Ny,Nu,0,tvropt);
end

%% Save
save(mfilename);