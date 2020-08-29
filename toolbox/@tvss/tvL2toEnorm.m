function [g,d,info] = tvL2toEnorm(G,Opt)
%% tvL2toEnorm
% This function computes the L2-to-Euclidean Norm exactly using the fact
% that it is a square root of the maximum eigenvalue of the matrix
% C(T)*P(T)*C(T)^T where P(T) is a controllability gramian evaluated at
% final time T.
%
% It uses the corresponding eigenvector to simulate the costate dynamics
% backwards in time to construct the worst-case disturbance that achives
% the specified gain.
%
% Initial conditions for plant must be zero. We will not get a clear
% expression of L2toE norm if initial conditions are nonzero.
%
% Inputs:
% G: LTV system defined on finite horizon
% Opt: [Optional] specfied using tvodeOptions

%% Input Processing
nin = nargin;
narginchk(1,2);
switch nin
    case 1
        Opt = tvodeOptions;
end

% Get Horizon
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
[~,~,CTf,DTf] = ssdata(tvsubs(G,Tf));

% Check Feedthrough
if any(DTf(:))
    error('The feedthrough matrix must be zero for L2toE gain to be well-defined.');
end

% Adjoint system
Ga = G';
Ga = Ga(:,[]);

% Boundary Condition
[A,B] = ssdata(G);
F = zeros(size(A,1));

%% Compute L2toE Gain & Construct worst-case disturbance
t0      = tic;
Wc      = cdle(A,B,F,[T0 Tf],Opt);
WcTf    = tvsubs(Wc,Tf);
M       = CTf*WcTf*CTf';
[V,D]   = eig(M);
[m,id]  = max(diag(D));
g       = sqrt(m);
V1      = V(:,id);
lamTf   = CTf'*V1;
sig     = tvmat(zeros(0,1),[T0,Tf]);
dwc     = tvlsim(Ga,sig,[Tf,T0],lamTf,Opt);
dwcnorm = tvnorm(dwc);
d       = dwc/dwcnorm;
tTotal  = toc(t0);

%% Process final output
info = [];
info.TotalTime = tTotal;
info.Wc = Wc;
end