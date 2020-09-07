function d = wcdist(A,B,gLow,PLow,COSTfh,T0,Nd,Opt)
%% Construct Worst-Case Disturbance Input From Incomplete RDE Solutions
% The construction is based on the two-point boundary value problem
% (TPBVP) related to the (Q,S,R,F) cost.  This particular implementation
% using the transformed Hamiltonian dynamics was developed by A. Iannelli
% and P. Seiler in September 2017.

% OdeSolver
OdeSolver = str2func(Opt.OdeSolver);
OdeOpt    = Opt.OdeOptions;

% Compute largest eigenvalue of PLow
T = PLow.Time;
[evecP,evalP] = eig( tvsubs(PLow,T(1)) );
[emax,idx] = max( diag(evalP) );
vmax = evecP(:,idx);

%norm( tvsubs(PLow,T(1))*vmax - emax*vmax )  %norm(emax*vmax)]
% Evaluate State/Cost Matrices on Same Time Grid as P
[~,R,S,~] = COSTfh(gLow);
[A,B,R,S] = evalt(A,B,tvmat(R),tvmat(S),T);

% Simulate Transformed Hamiltonian system
if isempty(S)
    M = R\(PLow*B)';
else
    M = R\(PLow*B+S)';
end
H11 = A-B*M;
odefh = @(t,E) tvsubs(H11,t)*E;
E0 = vmax/emax;
[tE,E] = OdeSolver( odefh,T,E0,OdeOpt );
E = tvmat(E,tE);

% Construct disturbance
% Note: This calculation uses a non-convergent Riccati solution and
% hence the disturbance starts at PLow.Time(1) > T0.
d = -M*E;

% Append zero input to disturbance on the window [PLow.Time(1) T0)
dTime = d.Time;
dData = d.Data;

dTime = [T0; T0+0.999*(dTime(1)-T0); dTime];
dData = cat(3,zeros(Nd,1,2), dData);
d = tvmat(dData,dTime);

% Normalize worst-case disturbance to have norm 1
d = d/tvnorm(d);
end