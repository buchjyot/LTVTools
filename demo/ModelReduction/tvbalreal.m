function [Gb,hSV,Ts,Tsi,info] = tvbalreal(G,Opt)
%% TVBALREAL Finite Horizon Gramian-based balancing of time-varying state-space realizations.
%
% [SYSB,G] = tvbalreal(SYS) computes a balanced state-space realization for
% the the linear time-varying system SYS. SYSB is an equivalent
% realization for which the controllability and observability Gramians are
% equal and diagonal, their diagonal entries forming the vector G of Hankel
% singular values. Small entries in G indicate states that can be removed
% to simplify the model (use MODRED to reduce the model order).
%
% [SYSB,G,T,Ti] = balreal(SYS,...) also returns the balancing state
% transformation x_b = T*x used to transform SYS into SYSB, as well as
% the inverse transformation x = Ti*x_b.

% Verify System is indeed finite horizon
[T0,Tf] = getHorizon(G); %#ok<ASGLU>
nin = nargin;
if isequal(nin,1)
    Opt = tvodeOptions('OdeSolver','ode23s');
end

% Compute finite horizon Gramians
[Wc] = tvgram(G,'c',Opt);
[Wo] = tvgram(G,'o',Opt);

% Unified Time Grid
Tgrid = union(union(Wc.Time,Wo.Time),G.Time);
[G,Wc,Wo] = evalt(G,Wc,Wo,Tgrid);

% Numerically Perturb boundry conditions
Nx = order(G);
Wc.Data(:,:,1) = diag(1e-10*ones(Nx,1));
Wo.Data(:,:,end) = diag(1e-10*ones(Nx,1));

% Compute Cholesky Factors
WcC = sqrtm(Wc);
WoC = sqrtm(Wo);
[U,S,V] = svd(WcC'*WoC);

% computing balanced transformations
hSV = diag(S);
sgi = 1./sqrt(hSV);
dsgi = diag(sgi);
Ts = dsgi*(U'*WoC);     % efficient diag(sgi)*u'*Ro
Tsi = (WcC'*V)*dsgi;  % efficient Rr'*v*diag(sgi)

% Balanced Realization
Ts.InterpolationMethod = 'Spline';
Tdot = tvdiff(Ts);
Ts.InterpolationMethod = 'Linear';
Tdot.InterpolationMethod = 'Linear';

[A,B,C,D] = ssdata(G);
Ab = (Ts*A + Tdot)*Tsi;
Bb = Ts*B;
Cb = C*Tsi;
Db = D;
Gb = tvss(Ab,Bb,Cb,Db);

% Process Output
info.Wc = Wc;
info.Wo = Wo;
info.Tdot = Tdot;
info.Tinv = Tsi;
end