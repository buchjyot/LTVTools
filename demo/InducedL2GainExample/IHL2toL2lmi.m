function [g,X11,P] = IHL2toL2lmi(G,v,p)
% This function computes the (infinite-horizon) induced L2 gain
% for an uncertain system Fu(G,Delta) using the IQC LMI condition.
% **The function assumes G is LTI and Delta is a SISO, unit norm bounded, 
% LTI uncertainty.
% 
% Inputs
%   G - Nominal LTI system
%   (v,p) - The IQC multiplier is blkdiag(Psiv'*X11*Psiv,-Psiv'*X11*Psiv)
%           where Psiv:=[1; 1/(s-p); ...; 1/(s-p)^v] and X11>=0.
% 
% Outputs
%   g - Gain Upper bound
%   X11 - IQC matrix
%   P - Storage function matrix

%% Form Extended System
[AA,BB,CC,DD,Nz,~,Nd] = ExtSystem(G,v,p);
n = size(AA,1);
CC1 = CC(1:Nz,1:n);
CC2 = CC(Nz+1:end,1:n);
DD1 = DD(1:Nz,:);
DD2 = DD(Nz+1:end,:);

%% Solve Dissipation Ineq LMI with IQC
% LMILab Implementation
setlmis([]);

% Create LMI variables
[X11,~,X11Dec] = lmivar(1,[v+1 1]);
X = lmivar(3,blkdiag(X11Dec,-X11Dec) );
P = lmivar(1,[n,1]);
[gsq,ndec] = lmivar(1,[1 1]);

% Define LMIs
lmiterm([1 1 1 X],[CC1 DD1]',[CC1 DD1]);
lmiterm([1 1 1 0],[CC2 DD2]'*[CC2 DD2]);
lmiterm([-1 1 1 gsq],[zeros(n+1,Nd); eye(Nd)],[zeros(n+1,Nd); eye(Nd)]');
lmiterm([1 1 1 P],[eye(n); zeros(Nd+1,n)],[AA BB],'s');

lmiterm([-2 1 1 P],1,1);
lmiterm([-3 1 1 X11],1,1);

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(end) = 1;

% Solve LMI
lmisys = getlmis;
opts = zeros(5,1);
%opts(1) = 1e-4;% Relative Accuracy
opts(2) = 200; % Max # of iters
opts(5) = 1;   % Toggle display
[copt,xopt] = mincx(lmisys,cobj,opts);
g = sqrt(copt);
X11 = dec2mat(lmisys,xopt,X11);
P = dec2mat(lmisys,xopt,P);
