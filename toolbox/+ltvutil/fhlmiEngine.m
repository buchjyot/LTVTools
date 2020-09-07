function [g,info] = fhlmiEngine(G,Delta,Ps,Psdot,Pm,Pmdot,Opt,NE)
% This function computes the finite-horizon induced L2 to Euclidean gain
% for an uncertain system Fu(G,Delta) using the IQC LMI condition.
% **The function assumes G is LTV and Delta is a SISO, unit norm bounded,
% LTI uncertainty.
%
% Inputs
%   G - Nominal LTV system (as a TVSS)
%   BlkStruct - Uncertainty Information
%   t - Nt-by-1 vector of time points to enforce the LMI
%   (Ps,Psdot) - Ns-by-1 vector of scalar basis functions and their
%        time derivatives (as TVMAT objects)
%   (Pm,Pmdot) - N-by-N matrix basis basis function and its time
%        derivative (as TVMAT objects)
%   NE - Number of outputs penalized in Euclidean sense
%
%    NOTE - All data (G,Ps,Psdot,Pm,Pmdot) must have the same time
%    grid t.  The LMI is enforced on this time grid t.
%
% Outputs
%   g - Gain Upper bound
%   X11 - IQC matrix
%   (P,Pdot) - Storage function matrix and its derivative.  P is returned
%         as an n-by-n-by-Nt array.

%% Rout the call to LOCAL function
flag = ltvutil.iqcEngine(Delta);
switch flag
    case 0
        % Enfore M11 > 0
        LMIFCN = @LOCALfhlmiDefault;
        
    case 1
        % Enfore M11(t) > 0
        LMIFCN = @LOCALfhlmiTVM11;
end

% Set NE to 0 if nargin < 8
if nargin < 8
    NE = 0;
end
[g,info] = LMIFCN(G,Delta,Ps,Psdot,Pm,Pmdot,Opt,NE);
end

function [g,info] = LOCALfhlmiDefault(G,Delta,Ps,Psdot,Pm,Pmdot,Opt,NE)
%% Input Processing
% XXX - Current allows only a single matrix basis function
% XXX - Verify that all objects have same time grid.
if isempty(Pm)
    useMatrixBasis = false;
else
    useMatrixBasis = true;
end
t = G.Time;
Ps = Ps(:);
Ns = size(Ps,1);
Nt = numel(t);

%% Read Uncertainty Information
[Nv,Nw] = size(Delta);
[~,Psi] = ltvutil.iqcEngine(Delta);
v = order(Psi)/2;

%% Form Extended System
[AAe,BBe,CCe,DDe,CCeT,Nz,~,Nd] = ltvutil.extsystem(G,Psi,Nv,Nw,NE);
n = size(AAe,1);
CC1e = CCe(1:Nz,1:n);
CC2e = CCe(Nz+1:end,1:n);
DD1e = DDe(1:Nz,:);
DD2e = DDe(Nz+1:end,:);

%% Solve Dissipation Ineq LMI with IQC
% LMILab Implementation
setlmis([]);

% Create LMI variables
[X11,~,X11Dec] = lmivar(1,[v+1 1]);
X = lmivar(3,blkdiag(X11Dec,-X11Dec) );
Y = zeros(Ns,1);
for k=1:Ns
    Y(k) = lmivar(1,[n 1]);
end
if useMatrixBasis
    lamP = lmivar(1,[1 1]);
end
[gsq,ndec] = lmivar(1,[1 1]);

% Define LMIs at each grid point
for i=1:Nt
    % Evaluate state matrices
    % XXX - This is not using any TVMAT methods, i.e. it assumes the
    % .Data subsref returns data in a Nr-by-Nc-by-Nt format.
    AA = AAe.Data(:,:,i);
    BB = BBe.Data(:,:,i);
    CC1 = CC1e.Data(:,:,i);
    CC2 = CC2e.Data(:,:,i);
    DD1 = DD1e.Data(:,:,i);
    DD2 = DD2e.Data(:,:,i);
    
    % Define LMIs:
    lmiterm([i 1 1 X],[CC1 DD1]',[CC1 DD1]);
    lmiterm([i 1 1 0],[CC2 DD2]'*[CC2 DD2]);
    lmiterm([-i 1 1 gsq],[zeros(n+1,Nd); eye(Nd)],[zeros(n+1,Nd); eye(Nd)]');
    Pspi = Ps.Data(:,:,i);
    Pspdoti = Psdot.Data(:,:,i);
    for j=1:Ns
        lmiterm([i 1 1 Y(j)],[Pspi(j)*eye(n); zeros(Nd+1,n)],...
            [AA BB],'s');
        lmiterm([i 1 1 Y(j)],0.5*[Pspdoti(j)*eye(n);zeros(Nd+1,n)],...
            [eye(n);zeros(Nd+1,n)]','s');
    end
    if useMatrixBasis
        Pmi = Pm.Data(:,:,i);
        Pmdoti = Pmdot.Data(:,:,i);
        lmiterm([i 1 1 lamP],[Pmi; zeros(Nd+1,n)],...
            [AA BB],'s');
        lmiterm([i 1 1 lamP],0.5*[Pmdoti;zeros(Nd+1,n)],...
            [eye(n);zeros(Nd+1,n)]','s');
    end
    
    % Add Terminal Penalties
    if i==Nt
        for j=1:Ns
            lmiterm([-(Nt+1) 1 1 Y(j)],0.5*Pspi(j)*eye(n),eye(n),'s');
        end
        if useMatrixBasis
            lmiterm([-(Nt+1) 1 1 lamP],0.5*Pmi,eye(n),'s');
        end
        if (NE~=0)
            lmiterm([(Nt+1) 1 1 0],CCeT'*CCeT);
        end
    end
end
lmiterm([-(Nt+2) 1 1 X11],1,1);

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(end) = 1;

% Solve LMI
lmisys = getlmis;
[copt,xopt] = mincx(lmisys,cobj,Opt.LmiSolverOptions);
if isempty(copt) && isempty(xopt)
    error('MINCX did not converge to a solution.');
end

gsq = copt;
X11 = dec2mat(lmisys,xopt,X11);
Yv = zeros(n,n,Ns);
for i=1:Ns
    Yv(:,:,i) = dec2mat(lmisys,xopt,Y(i));
end
Y = Yv;

if useMatrixBasis
    lamP = dec2mat(lmisys,xopt,lamP);
else
    lamP = [];
end

%% Store Outputs
g = sqrt(gsq);
X11 = full(X11);

P = Ps(1)*Y(:,:,1);
Pdot = Psdot(1)*Y(:,:,1);
for j=2:Ns
    P = P + Ps(j)*Y(:,:,j);
    Pdot = Pdot + Psdot(j)*Y(:,:,j);
end
if useMatrixBasis
    % XXX
    Pm.InterpolationMethod = 'Spline';
    Pmdot.InterpolationMethod = 'Spline';
    
    P = P + lamP*Pm;
    Pdot = Pdot + lamP*Pmdot;
end

info = struct('X11',X11,'P',P,'Pdot',Pdot,'Y',Y,'lamP',lamP);
end

function [g,info] = LOCALfhlmiTVM11(G,Delta,Ps,Psdot,Pm,Pmdot,Opt,NE)
% XXX - Assumes Delta is a memoryless, TV function unit norm-bounded.
% Code below is modified to allow IQC matrix to be time-varying.
% XXX - v=0 is required here.....
%
% This function computes the finite-horizon induced L2 to Euclidean gain
% for an uncertain system Fu(G,Delta) using the IQC LMI condition.
% **The function assumes G is LTV and Delta is a SISO, unit norm bounded,
% LTI uncertainty.
%
% Inputs
%   G - Nominal LTV system (as a TVSS)
%   (v,p) - The IQC multiplier is blkdiag(Psiv'*X11*Psiv,-Psiv'*X11*Psiv)
%           where Psiv:=[1; 1/(s-p); ...; 1/(s-p)^v] and X11>=0.
%   t - Nt-by-1 vector of time points to enforce the LMI
%   (Ps,Psdot) - Ns-by-1 vector of scalar basis functions and their
%        time derivatives (as TVMAT objects)
%   (Pm,Pmdot) - N-by-N matrix basis basis function and its time
%        derivative (as TVMAT objects)
%   NE - Number of outputs penalized in Euclidean sense
%
%    NOTE - All data (G,Ps,Psdot,Pm,Pmdot) must have the same time
%    grid t.  The LMI is enforced on this time grid t.
%
% Outputs
%   g - Gain Upper bound
%   X11 - IQC matrix (XXX - TVMAT)
%   (P,Pdot) - Storage function matrix and its derivative.  P is returned
%         as an n-by-n-by-Nt array.

%% Input Processing
% XXX - Current allows only a single matrix basis function
% XXX - Verify that all objects have same time grid.
if isempty(Pm)
    useMatrixBasis = false;
else
    useMatrixBasis = true;
end
t = G.Time;
Ps = Ps(:);
Ns = size(Ps,1);
Nt = numel(t);

%% Read Uncertainty Information
[Nv,Nw] = size(Delta);
[~,Psi] = ltvutil.iqcEngine(Delta);
v = order(Psi)/2;
if v~=0
    error('This function requires filter order to be 0.');
end

%% Form Extended System
[AAe,BBe,CCe,DDe,CCeT,Nz,~,Nd] = ltvutil.extsystem(G,Psi,Nv,Nw,NE);
n = size(AAe,1);
CC1e = CCe(1:Nz,1:n);
CC2e = CCe(Nz+1:end,1:n);
DD1e = DDe(1:Nz,:);
DD2e = DDe(Nz+1:end,:);

%% Solve Dissipation Ineq LMI with IQC
% LMILab Implementation
setlmis([]);

% Create LMI variables
X11 = zeros(Nt,1);
X = zeros(Nt,1);
for k=1:Nt
    [X11(k),~,X11Dec] = lmivar(1,[v+1 1]);
    X(k) = lmivar(3,blkdiag(X11Dec,-X11Dec) );
end

Y = zeros(Ns,1);
for k=1:Ns
    Y(k) = lmivar(1,[n 1]);
end
if useMatrixBasis
    lamP = lmivar(1,[1 1]);
end
[gsq,ndec] = lmivar(1,[1 1]);

% Define LMIs at each grid point
for i=1:Nt
    % Evaluate state matrices
    % XXX - This is not using any TVMAT methods, i.e. it assumes the
    % .Data subsref returns data in a Nr-by-Nc-by-Nt format.
    AA = AAe.Data(:,:,i);
    BB = BBe.Data(:,:,i);
    CC1 = CC1e.Data(:,:,i);
    CC2 = CC2e.Data(:,:,i);
    DD1 = DD1e.Data(:,:,i);
    DD2 = DD2e.Data(:,:,i);
    
    % Define LMIs:
    lmiterm([i 1 1 X(i)],[CC1 DD1]',[CC1 DD1]);
    lmiterm([i 1 1 0],[CC2 DD2]'*[CC2 DD2]);
    lmiterm([-i 1 1 gsq],[zeros(n+1,Nd); eye(Nd)],[zeros(n+1,Nd); eye(Nd)]');
    Pspi = Ps.Data(:,:,i);
    Pspdoti = Psdot.Data(:,:,i);
    for j=1:Ns
        lmiterm([i 1 1 Y(j)],[Pspi(j)*eye(n); zeros(Nd+1,n)],[AA BB],'s');
        lmiterm([i 1 1 Y(j)],0.5*[Pspdoti(j)*eye(n);zeros(Nd+1,n)],...
            [eye(n);zeros(Nd+1,n)]','s');
    end
    if useMatrixBasis
        Pmi = Pm.Data(:,:,i);
        Pmdoti = Pmdot.Data(:,:,i);
        lmiterm([i 1 1 lamP],[Pmi; zeros(Nd+1,n)],[AA BB],'s');
        lmiterm([i 1 1 lamP],0.5*[Pmdoti;zeros(Nd+1,n)],...
            [eye(n);zeros(Nd+1,n)]','s');
    end
    
    if i==Nt
        for j=1:Ns
            lmiterm([-(Nt+1) 1 1 Y(j)],0.5*Pspi(j)*eye(n),eye(n),'s');
        end
        if useMatrixBasis
            lmiterm([-(Nt+1) 1 1 lamP],0.5*Pmi,eye(n),'s');
        end
        if (NE~=0)
            lmiterm([(Nt+1) 1 1 0],CCeT'*CCeT);
        end
    end
end

for k=1:Nt
    lmiterm([-(Nt+1+k) 1 1 X11(k)],1,1);
end

% Create objective function: min gam
cobj = zeros(ndec,1);
cobj(end) = 1;

% Solve LMI
lmisys = getlmis;
[copt,xopt] = mincx(lmisys,cobj,Opt.LmiSolverOptions);
if isempty(copt) && isempty(xopt)
    error('MINCX did not converge to a solution.');
end

gsq = copt;
X11Data = zeros(Nt,1);
for i=1:Nt
    X11Data(i) = dec2mat(lmisys,xopt,X11(i));
end
X11 = tvmat(X11Data,t);

Yv = zeros(n,n,Ns);
for i=1:Ns
    Yv(:,:,i) = dec2mat(lmisys,xopt,Y(i));
end
Y = Yv;

if useMatrixBasis
    lamP = dec2mat(lmisys,xopt,lamP);
else
    lamP = [];
end

%% Store Outputs
g = sqrt(gsq);
% X11 = full(X11);

P = Ps(1)*Y(:,:,1);
Pdot = Psdot(1)*Y(:,:,1);
for j=2:Ns
    P = P + Ps(j)*Y(:,:,j);
    Pdot = Pdot + Psdot(j)*Y(:,:,j);
end
if useMatrixBasis
    % XXX
    Pm.InterpolationMethod = 'Spline';
    Pmdot.InterpolationMethod = 'Spline';
    
    P = P + lamP*Pm;
    Pdot = Pdot + lamP*Pmdot;
end
info = struct('X11',X11,'P',P,'Pdot',Pdot,'Y',Y,'lamP',lamP);
end