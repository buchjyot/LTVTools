function [g,X11,P,Pdot,Y,lamP] = FHL2toL2lmi(G,v,p,Ps,Psdot,Pm,Pmdot)
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
%
%    NOTE - All data (G,Ps,Psdot,Pm,Pmdot) must have the same time
%    grid t.  The LMI is enforced on this time grid t.
%
% Outputs
%   g - Gain Upper bound
%   X11 - IQC matrix
%   (P,Pdot) - Storage function matrix and its derivative.  P is returned
%         as an n-by-n-by-Nt array.

%% Input Processing
% XXX - Current allows only a single matrix basis function
% XXX - Verify that all objects have same time grid.
% XXX - G must have zero feedthrough from d to e for L2 to E. Check this.
if nargin==5 || isempty(Pm)
    useMatrixBasis = false;
else
    useMatrixBasis = true;
end
t = G.Time;
Ps = Ps(:);
Ns = size(Ps,1);
Nt = numel(t);

%% Form Extended System
[AAe,BBe,CCe,DDe,Nz,~,Nd] = ExtSystem(G,v,p);
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
    ti = t(i);
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
    Psi = Ps.Data(:,:,i);
    Psdoti = Psdot.Data(:,:,i);
    for j=1:Ns
        lmiterm([i 1 1 Y(j)],[Psi(j)*eye(n); zeros(Nd+1,n)],...
            [AA BB],'s');
        lmiterm([i 1 1 Y(j)],0.5*[Psdoti(j)*eye(n);zeros(Nd+1,n)],...
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
    
    if i==Nt
        for j=1:Ns
            lmiterm([-(Nt+1) 1 1 Y(j)],0.5*Psi(j)*eye(n),eye(n),'s');
        end
        if useMatrixBasis
            lmiterm([-(Nt+1) 1 1 lamP],0.5*Pmi,eye(n),'s');
        end
    end
end
lmiterm([-(Nt+2) 1 1 X11],1,1);

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

