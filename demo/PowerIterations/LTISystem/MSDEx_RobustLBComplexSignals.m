%% Mass-Spring-Damper System
% m := mass
% b := damping coefficient
% k := spring constant
msdEx1 = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
P  = msdEx1(1,2*0.7,1);

% Transform plant to have all the states (Nx) as outputs
Nx = order(P);
P = ss(P.A,P.B,eye(Nx),0);
[Ne,Nd] = size(P);

%% Horizon
T0   = 0;
Tall = 3;%[0.1,1,2,3,4,5,7,10,20,50,70];
NT   = length(Tall);

%% Analysis Interconnection
% Input Scalling
WdScl = 0.1;

% Output Scalling
WexScl = 1;

% Uncertainty Scalling
Nv = 1;
Nw = 1;
DelNorm = 0.6;
DelInScl = sqrt(DelNorm);
DelOutScl = sqrt(DelNorm);

% Promote to state-space objects
Wd = ss(WdScl*eye(Nd));
We = ss(WexScl*eye(Ne));
WDelIn = ss(DelInScl*eye(Nv));
WDelOut = ss(DelOutScl*eye(Nw));

% System Interconnection
systemnames = 'P Wd We WDelIn WDelOut';
inputvar = '[w; d]';
outputvar = '[WDelIn; We]';
input_to_P = '[WDelOut+Wd]';
input_to_Wd = '[d]';
input_to_We = '[P]';
input_to_WDelIn = '[Wd]';
input_to_WDelOut = '[w]';
cleanupsysic = 'yes';
Gunc = sysic;

% Euclidean Penalties
NE = 0;

%% Options
OdeSolver = 'ode23s';
pOpt  = poweritOptions('Display','on','OdeSolver',OdeSolver,'StoreAllIter',true);
pSpec = poweritSignalSpec('Nv',1,'Nw',1,'NE',NE,'InitialInput','randnc');

% seed
rng(0);

%% Power Iterations Lower Bound
ComputeWcGainLowerBounds = true;
if ComputeWcGainLowerBounds
    wcgLB = zeros(NT,1);
    pInfo = cell(NT,1);
    wcSig = cell(NT,1);
    for i = 1:NT
        [wcgLB(i),wcSig{i},pInfo{i}] = powerit(Gunc,[T0,Tall(i)],pSpec,pOpt);
        fprintf(newline);
    end
end

%% Save Data
save(mfilename);
return;

%% Analysis for i = 1
Gt = tvss(Gunc,[0,Tall(1)]);

% Verify worst-case gain for computed input
Uwc = wcSig{1};

% wcgSim should match with wcgLB
[wcgSim,Ysim,Xsim] = tvlsimwcgain(Gt,Uwc,1,1,NE);

% Read inputs
w = Uwc(1);
u = Uwc(2);

% Read worstcase outputs
v = Ysim(1);
yL2 = Ysim(2:3);

% Construct Uncertainty