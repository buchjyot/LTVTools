function [g,d,info] = tvnorm2(G,gTry,NE,Opt)
%% TVNORM2 Checks if specified gTry is an upper bound or lower bound
%
% This is a special function using which you can check if specified gTry is
% an upper bound or lower bound for the specified LTV system. Here, the
% cost is specified through number of Euclidean outputs NE.
%
% NOTE:
% Let G be a TVSS system with NY by NU IO dimentions
% NE is an optional argument with default set to 0.
% NE = NY means 'L2toE' Norm
% 0 < NE < NY then we have 'L2toE' and 'L2toL2' mixed cost terms.
%
% Syntex:
% >> [g,d,info] = tvnorm2(G,gTry,NE,Opt);
%
% If gTry is a lower bound then output g = [gTry NaN] and we return the
% worst-case disturbance constructed using whatever RDE solution we have.
%
% If gTry is an upper bound then output g = [NaN gTry] and d = [].
%
% Requires at least 2 input arguments, 3rd input argument is NE and 4th
% argument is Opt, which are both optional. NE = 0 is a default. Opt can be
% specified through tvodeOptions, we do not require tvnormOptions, because
% no bisection is performed in this function.

% Input Processing
narginchk(2,4);
nin = nargin;
nout = nargout;
switch nin
    case 2
        NE = 0;
        Opt = tvodeOptions;
    case 3
        Opt = tvodeOptions;
end

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
NY = size(G,1);

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY.');
end

%% Build Cost Function Matrices
[A,B] = ssdata(G);
[Nx,Nd] = size(B);

% If Nx == 0 then error out
if isequal(Nx,0)
    error('The state dimention must be nonzero.');
end

% We do not support discriptor systems currently
E = [];

% XXX (JB) The following should work when 1:NY-NE is empty by some size
% >> G(1:NY-NE,:)

% Extract L2 part
[~,~,CL2,DL2] = ssdata(tvss(G.Data(1:NY-NE,:),G.Time,G.InterpolationMethod));

% Extract Euclidean part
[~,~,CE,DE] = ssdata(tvss(G.Data(NY-NE+1:NY,:),G.Time,G.InterpolationMethod));
[CTf,DTf] = tvsubs(CE,DE,Tf);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm');
end

% Generic Cost Matrices
% XXX This set-up requires R(t) to be inverted at each time step.
% However, no inversion is actually required since we can set
% R=I and B'=B/gamma.  This should speed up the integration.
R0 = DL2'*DL2;
R1 = eye(Nd);
R = R0-(gTry^2)*R1;
Q = CL2'*CL2;
S = CL2'*DL2;
F = CTf'*CTf;

% Place holders for time recording
tCDRE = []; %#ok<NASGU>
tWCDistConstruction = [];

%% Solve CDRE
% Record time required for CDRE integration
t0 = tic;
if nout==1
    P = cdre(A,B,Q,R,S,E,F,[Tf T0],Opt);
    Pdot = []; sol = [];
else
    [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],Opt);
end
tCDRE = toc(t0);

% Check convergence of P
converged_flag = ~(P.Time(1)>T0);
if converged_flag
    % P converged (gTry is an upper bound)
    g = [NaN gTry];
else
    % P did not converge (gTry is a lower bound)
    g = [gTry NaN];
end

%% Construct Worst-Case Input
% The construction is based on the two-point boundary value problem
% (TPBVP) related to the (Q,S,R,F) cost.  This particular implementation
% using the transformed Hamiltonian dynamics was developed by A. Iannelli
% and P. Seiler in September 2017.
d = [];
if nout>=2 && ~converged_flag
    % Start timing
    t1 = tic;
    
    % Compute largest eigenvalue of P
    T = P.Time;
    [evecP,evalP] = eig( tvsubs(P,T(1)) );
    [emax,idx] = max( diag(evalP) );
    vmax = evecP(:,idx);
    
    %norm( tvsubs(P,T(1))*vmax - emax*vmax )  %norm(emax*vmax)]
    % Evaluate State/Cost Matrices on Same Time Grid as P
    [~,R,S,~] = COSTfh(gLow);
    [A,B,R,S] = evalt(A,B,tvmat(R),tvmat(S),T);
    
    % Simulate Transformed Hamiltonian system
    if isempty(S)
        M = R\(P*B)';
    else
        M = R\(P*B+S)';
    end
    H11 = A-B*M;
    odefh = @(t,E) tvsubs(H11,t)*E;
    E0 = vmax/emax;
    [tE,E] = ode45( odefh,T,E0);
    E = tvmat(E,tE);
    
    % Construct disturbance
    % Note: This calculation uses a non-convergent Riccati solution and
    % hence the disturbance starts at P.Time(1) > T0.
    d = -M*E;
    
    % Append zero input to disturbance on the window [P.Time(1) T0)
    dTime = d.Time;
    dData = d.Data;
    
    dTime = [T0; T0+0.999*(dTime(1)-T0); dTime];
    dData = cat(3,zeros(Nd,1,2), dData);
    d = tvmat(dData,dTime);
    
    % Normalize worst-case disturbance to have norm 1
    d = d/tvnorm(d);
    
    % Record time
    tWCDistConstruction = toc(t1);
end

%% Store Final Result
if nout>=3
    info = [];
    info.P = P;
    info.Pdot = Pdot;
    info.sol  = sol;
    info.Gain = g;
    info.IntegrationTime = tCDRE;
    info.ConstructWCDistTime = tWCDistConstruction;
end
end