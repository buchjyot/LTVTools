function [g,d,info] = tvnorm1(G,varargin)
%% TVNORM Computes the time varying system norm
%
% This function combines the power iteration and bisection method.
%
% [g,d,info] = tvnorm1(G);
% Computes default 'L2toL2' norm with default options.
%
% [g,d,info] = tvnorm1(G,Opt);
% Computes default 'L2toL2' with user specified options.
%
% [g,d,info] = tvnorm1(G,NE);
% Computes tvnorm based on following cases with default options.
%
% [g,d,info] = tvnorm1(G,NE,Opt);
% Computes tvnorm based on following cases with user specified Options.
%
% NOTE:
% Let G be a TVSS system with NY by NU IO dimentions
% NE is an optional argument with default set to 0.
% NE = NY means 'L2toE' Norm
% 0 < NE < NY then we have 'L2toE' and 'L2toL2' mixed cost terms.
%
% Output Arguments are worst-case gain (g) induced by L2 disturbance (d)
% and relevant Riccati solution info

%% Input Processing

% Check # of inputs
narginchk(1,3);
nin = nargin;
nout = nargout;

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
NY = size(G,1);

% Obtain default Opt
switch nin
    case 1
        Opt = tvnormOptions;
        NE = 0;
    case 2
        if isa(varargin{1},'tvnormOptions')
            Opt = varargin{1};
            NE = 0;
        else
            Opt = tvnormOptions;
            NE = varargin{1};
        end
    case 3
        NE = varargin{1};
        Opt = varargin{2};
end

% DispFlag
DispFlag = isequal(Opt.Display,'on');

% Update cdreOpt based on tvnormoptions
odeOpt = tvodeOptions('OdeOptions',Opt.OdeOptions,'OdeSolver',Opt.OdeSolver);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% Counter for RDE bisections
RDEcnt = 0;

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY.');
end

%% Build Cost Function Matrices
[A,B]   = ssdata(G);
[Nx,Nd] = size(B);
NL2     = NY-NE;

% If Nx == 0 then error out
if isequal(Nx,0)
    error('The state dimention must be nonzero.');
end

% We do not support discriptor systems currently
E = [];

% XXX (JB) The following should work when 1:NY-NE is empty by some size
% >> G(1:NY-NE,:)

% Extract L2 part
[~,~,CL2,DL2] = ssdata(tvss(G.Data(1:NL2,:),G.Time,G.InterpolationMethod));

% Extract Euclidean part
[~,~,CE,DE] = ssdata(tvss(G.Data(NL2+1:NY,:),G.Time,G.InterpolationMethod));
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
Rc = @(GAM) R0-(GAM^2)*R1;
Qc = CL2'*CL2;
Sc = CL2'*DL2;
Fc = CTf'*CTf;
COSTfh = @(GAM) deal(Qc,Rc(GAM),Sc,Fc);

%% Initialization
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];

AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;

% poweritOptions
Opt = tvpoweritOptions('StopTol',5e-3,'Display',Opt.Display);

%% Lower Bound Phase
% Power iterations
if DispFlag
    fprintf(' ### Starting Power Iterations:\n');
end
t1 = tic;
GTime = G.Time;
Nt = length(GTime);
d1 = tvmat(randn(Nd,1,Nt),GTime);
d1 = d1/tvnorm(d1);
[gLow,dwc] = LOCALPowerit(G,Fc,Tf,NE,NL2,NY,Nx,d1,Opt);
if DispFlag
    fprintf(' Lower Bound = %4.3f\n',gLow);
end

%% Upper Bound Phase
% Attempt to determine a finite upper bound on performance.
gFac = 1.1; % Keep increasing 10% until you find an upper bound
gUpp = (2*(gLow+AbsTol+(gLow*0.05/100))/(1-RelTol)) - gLow;  % XXX Better choice?

cnt = 0;
cntMax = 8;
haveUpper = false;
while ~haveUpper && cnt<cntMax
    % Pick gamma
    cnt = cnt+1;
    gTry = gUpp;
    
    % Solve LTV Riccati Equation
    [Q,R,S,F] = COSTfh(gTry);
    if nout==1
        P = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
        Pdot = []; sol = [];
    else
        [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
    end
    
    % Check convergence of P
    if P.Time(1)>T0
        % P did not converge
        gUpp = gFac*gUpp;
        gLow = gTry;
        PLow = P; PdotLow = Pdot; solLow = sol;
    else
        % P converged
        haveUpper = true;
        gUpp = gTry;
        PUpp = P; PdotUpp = Pdot; solUpp = sol;
    end
end

if ~haveUpper
    % Could not find a finite upper bound
    if DispFlag
        fprintf([' Could not find a finite upper bound.' ...
            ' Infeasible at gTry = %4.1f\n'],gUpp);
    end
    gUpp = inf;
end

if DispFlag
    fprintf(' Upper Bound = %4.3f\n',gUpp);
end

%% Bisection Phase
if DispFlag
    fprintf(' ### Starting Bisection Phase:\n');
end

if isfinite(gUpp)
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        gTry = (gUpp+gLow)/2;
        
        % Solve LTV Riccati Equation
        [Q,R,S,F] = COSTfh(gTry);
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
        end
        RDEcnt = RDEcnt + 1;
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gLow = gTry;
            PLow = P; PdotLow = Pdot; solLow = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        else
            % P converged
            gUpp = gTry;
            PUpp = P; PdotUpp = Pdot; solUpp = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        end
    end
    if isequal(RDEcnt,0) && DispFlag
        fprintf(' No bisection required.\n');
    end
end

%% Construct Worst-Case Input
% The construction is based on the two-point boundary value problem
% (TPBVP) related to the (Q,S,R,F) cost.  This particular implementation
% using the transformed Hamiltonian dynamics was developed by A. Iannelli
% and P. Seiler in September 2017.
d = [];
if nout>=2
    if ~isempty(PLow)
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
    else
        % Set the disturbance that we get from power iteration
        d = dwc;
    end
end

%% Store Final Result
tTotal = toc(t1);
g = [gLow,gUpp];

info.Lower.Gain = gLow;
info.Lower.P    = PLow;
info.Lower.Pdot = PdotLow;
info.Lower.sol  = solLow;

info.Upper.Gain = gUpp;
info.Upper.P    = PUpp;
info.Upper.Pdot = PdotUpp;
info.Upper.sol  = solUpp;

info.TotalBisections = RDEcnt;
info.TotalTime       = tTotal;
end

function [glb,dwc] = LOCALPowerit(G,F,Tf,NE,NL2,Ny,Nx,d1,Opt)
%% LOCALPowerit
% Runs power iterations to estimate the best lower bound

% Adjoint System
Ga = G';

% Initial Condition
x0 = zeros(Nx,1);

% Options
MaxIter         = Opt.MaxIter;
StopTol         = Opt.StopTol;
OdeSolver       = Opt.OdeSolver;
OdeOptions      = Opt.OdeOptions;
StepSize        = Opt.StepSize;

% Options
Display         = isequal(Opt.Display,'on');
tvsopt1         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','ForwardInTime','StepSize',StepSize);
tvsopt2         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','BackwardInTime','StepSize',StepSize);

% Memory Allocation
U    = cell(MaxIter,1);
Y    = cell(MaxIter,1);
X    = cell(MaxIter,1);
W    = cell(MaxIter,1);
Perf = zeros(MaxIter,1);

% Display String & Performance
DisplayProgress = @(id,pf) fprintf(' Iter: %d, Performance: %.3f\n',id,pf);
if NE == 0
    % L2toL2
    EvaluatePerformance = @(Yi,Xi) tvnorm(Yi);
    lamTf               = @(Xi) zeros(Nx,1);
    getCostateInputs    = @(Yi,pf) Yi/pf;
elseif NE == Ny
    % L2toE
    EvaluatePerformance = @(Yi,Xi) sqrt(tvsubs(Xi,Tf)'*F*tvsubs(Xi,Tf));
    lamTf               = @(Xi) F*tvsubs(Xi,Tf);
    getCostateInputs    = @(Yi,pf) evalt(tvmat(zeros(NE,1)),Yi.Time);
else
    % Mixed Cost
    EvaluatePerformance = @(Yi,Xi) sqrt(tvsubs(Xi,Tf)'*F*tvsubs(Xi,Tf) + tvnorm(Yi(1:NL2))^2);
    lamTf               = @(Xi) F*tvsubs(Xi,Tf);
    getCostateInputs    = @(Yi,pf) [Yi(1:NL2)/tvnorm(Yi(1:NL2)); zeros(NE,1)];
end

% Initial Input
U{1} = d1;

% Power Iterations
for i = 1:MaxIter
    
    % State Equation
    [Y{i},X{i}] = tvlsim(G,U{i},x0,tvsopt1);
    
    % Evaluate Performance
    Perf(i) = EvaluatePerformance(Y{i},X{i});
    
    % Display
    if Display
        DisplayProgress(i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter-1)
        % Check if performance that we care is stationary
        if abs(Perf(i)-Perf(i-1)) <= StopTol*max(Perf([i i-1]))
            break;
        end
    end
    
    % Get Boundary Conditions
    lambdaTf = lamTf(X{i});
    
    % Get Adjoint System Forcing (Alignment Condition for Y)
    W{i} = getCostateInputs(Y{i},Perf(i));
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,W{i},lambdaTf,tvsopt2);
    
    % Alignment Condition for U
    U{i+1} = U{i+1}/tvnorm(U{i+1});
end

% Final Output
[mPerf,id] = max(Perf(1:i));
glb = mPerf;
dwc = U{id};
end