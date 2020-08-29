%% Load LTV model data
load('twoLinkRobot_BuildLTVModel.mat','G','DelNorm','T0','Tf');
load('twoLinkRobot_LQR','Klqr');

% Evaluate Klqr and G on the same time grid
[G,Klqr] = evalt(G,Klqr,union(G.Time,Klqr.Time));

% Uncertain CL System:  Tunc = Fu(Tnom,Delta)
systemnames = 'G Klqr';
inputvar = '[w; d(2)]';
outputvar = '[d(2)-Klqr(2); G(1:2)]';
input_to_G = '[d(1)-Klqr(1); d(2)-Klqr(2)+w]';
input_to_Klqr = '[G]';
cleanupsysic = 'yes';
Tnom = sysic;

% Uncertainity Scalling
OutScl = blkdiag(sqrt(DelNorm),eye(2));
InScl  = blkdiag(sqrt(DelNorm),eye(2));
Tnom = OutScl*Tnom*InScl;

% Power-iterations lower bound
pSpec = poweritSignalSpec('NE',2,'Nv',1,'Nw',1);
pOpt = poweritOptions('Display','on');
[glb,dwc,info] = powerit(Tnom,[T0,Tf],pSpec,pOpt);