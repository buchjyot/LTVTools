function [g,d,info] = tvL2toE(A,B,CTf,Nx,T0,Tf,nout,Opt)
%% tvL2toE
% This function computes the L2-to-Euclidean Norm exactly using the fact
% that it is a square root of the maximum eigenvalue of the matrix
% C(T)*P(T)*C(T)^T where P(T) is a controllability gramian evaluated at
% final time T.
%
% It uses the corresponding eigenvector to simulate the costate dynamics
% backwards in time to construct the worst-case disturbance that achives
% the specified gain.

%% Input Processing
% Adjoint system
Ga = tvss(-A',[],B',[]);

% Boundary Condition
F = zeros(Nx);

% Initialize
d = []; info = [];

%% Compute L2toE Gain & Construct worst-case disturbance
t0      = tic;
Wc      = cdle(A,B,F,[T0 Tf],Opt);
WcTf    = tvsubs(Wc,Tf);
M       = CTf*WcTf*CTf';
[V,D]   = eig(M);
[m,id]  = max(diag(D));
g       = sqrt(m);

% Construct disutbance only if nout > 1
if nout > 1
    V1      = V(:,id);
    lamTf   = CTf'*V1;
    dwc     = tvlsim(Ga,[],[Tf,T0],lamTf,Opt);
    dwcnorm = tvnorm(dwc);
    d       = dwc/dwcnorm;
end

% Store info
info.TotalTime = toc(t0);
info.Wsc = Wc;
info.Woc = M;
if nout > 1
    info.V1 = V1;
end
end