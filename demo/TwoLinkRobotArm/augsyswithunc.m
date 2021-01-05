function [sysAug,sysAdj] = augsyswithunc(Gunc,Delta,vbar,NE)
%% AugSysWithUnc Constructs Augmented System with Uncertainty Delta and Affine Term
% Given:
%   Uncertainy system Gunc
%   Uncertainty delta
%   Additional input vbar
%   Number of Euclidean outputs
%
% Output: 
%   Augmented LTV system and its adjoint 
%
% NOTE : The Augmented System has non-zero initial conditions, See
% sysAug.UserData for this non-zero initial condition.

nout = nargout;

% Verify G is finite horizon
ltvutil.verifyFH(Gunc)

% Get uncertainty sizes
[~,Nv] = size(Delta);

% Make sure size of vbar is of same size as Nv
Nvbar = size(vbar,1);
if ~isequal(Nvbar,Nv)
   error('Trim input vbar must be TVMAT of size Nv by 1.') 
end

% Perform Interconnection and Construct Augmented system
temp        = lft( Delta, Gunc.Data );
G           = tvss(temp, Gunc.Time);
[A,B,C,D]   = ssdata(G);
Nx          = order(G);
[NY,NU]     = size(G);
tVec        = G.Time;

% Separate B and D
B1 = B(:,1:NU-Nv);
B2 = B(:,NU-Nv+1:end);
D1 = D(:,1:NU-Nv);
D2 = D(:,NU-Nv+1:end);

% build Abar
Azeros = tvmat(zeros(1,Nx+1,numel(tVec)),tVec);
Abar = [A, B2*vbar; Azeros];

% build Bbar
Bzeros = tvmat(zeros(1,NU-Nv,numel(tVec)),tVec);
Bbar = [B1; Bzeros];

% build Cbar
Cbar = [C, D2*vbar];

% build Dbar
Dbar = D1;

% System and Adjoint
sysAug = tvss(Abar,Bbar,Cbar,Dbar);
sysAug.UserData = [zeros(Nx,1); ones(Nv,1)];
if nout>1
    NL2 = NY-NE;
    sysAdj = tvss(-Abar',-Cbar(1:NL2,:)',Bbar',Dbar(1:NL2,:)');
end