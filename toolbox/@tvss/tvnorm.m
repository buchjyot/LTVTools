function [g,d,info] = tvnorm(G,varargin)
%% TVNORM
% Computes the time varying system norm using the approach in:
%
% Jyot Buch, Murat Arcak, Peter Seiler, An Efficient Algorithm to Compute
% Norms for Finite Horizon, Linear Time-Varying System, Submitted to
% LCSS/ACC 2021.
%
% [g,d,info] = tvnorm(G);
% Computes default 'L2toL2' norm with default options.
%
% [g,d,info] = tvnorm(G,Opt);
% Computes default 'L2toL2' with user specified options.
%
% [g,d,info] = tvnorm(G,NE);
% Computes tvnorm based on following cases with default options.
%
% [g,d,info] = tvnorm(G,NE,Opt);
% Computes tvnorm based on following cases with user specified Options.
%
% NOTE:
% Let G be a TVSS system with NY by NU IO dimentions
% NE is an optional argument with default set to 0.
% NE = NY means 'L2toE' Norm
% 0 < NE < NY then we have 'L2toE' and 'L2toL2' mixed cost terms.

%% Input Processing

% Check # of inputs
narginchk(1,3);
nout = nargout;

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
[NY,NU] = size(G);

% ParseInputs
[NE,Opt] = ltvutil.parseInput(varargin,{'double','tvnormOptions'});

% Default
if isempty(NE)
    NE = 0;
end
if isempty(Opt)
    Opt = tvnormOptions;
end

% DispFlag
Display = isequal(Opt.Display,'on');

% Update cdreOpt based on tvnormoptions
odeOpt = tvodeOptions('OdeOptions',Opt.OdeOptions,'OdeSolver',Opt.OdeSolver);

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY.');
end

%% Build Cost Function Matrices
[A,B,C,D]   = ssdata(G);
[Nx,Nd]     = size(B);
NL2         = NY-NE;

% If Nx == 0 then error out
if isequal(Nx,0)
    error('The state dimention must be nonzero.');
end

% We do not support discriptor systems
E = [];

% Extract Euclidean part
[CE,DE] = tvsubs(C,D,Tf);
CTf = CE(NL2+1:NY,:);
DTf = DE(NL2+1:NY,:);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm');
end

% If NL2 == 0 then call utility function to compute L2toE
if isequal(NL2,0)
    [g,d,info] = ltvutil.tvL2toE(A,B,CTf,Nx,T0,Tf,nout,odeOpt);
    return;
end

% Extract L2 part
CL2 = C(1:NL2,:);
DL2 = D(1:NL2,:);

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

% Precompute Adjoint System
Ga = tvss(-A',-CL2',B',DL2');

%% Initialization
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];

AbsTol     = Opt.AbsTol;
RelTol     = Opt.RelTol;
gLow       = Opt.Bounds(1);
gUpp       = Opt.Bounds(2);
GTime      = G.Time;
Nt         = max(20,length(GTime));
x0         = zeros(Nx,1);
RDEcnt     = 0;
MaxPowIter = 50;
piterInfo  = cell(20,1);
d          = tvmat(randn(NU,1,Nt),linspace(T0,Tf,Nt));
d          = d/tvnorm(d);

%% Algorithm
% Start global timing
t0 = tic;

% Main loop
while (gUpp - gLow >= RelTol*gUpp + AbsTol)
    
    % Power iteration lower bound
    [gLow,d,piterInfo{RDEcnt+1}] = ltvutil.tvpowerit(G,Ga,d,CTf,T0,Tf,NL2,x0,MaxPowIter,RelTol/5,Display,odeOpt);
    
    % Try cleverly chosen bound
    gTry = (gLow + AbsTol)/(1 - RelTol);
    
    % Solve LTV Riccati Equation using gTry
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
        if Display
            fprintf(' gTry = %4.3f \t tConv = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                gTry,P.Time(1),gLow,gUpp);
        end
        
        % Construct Worst-Case disturbance
        d = ltvutil.wcdist(A,B,gLow,PLow,COSTfh,T0,Nd,Opt);
    else
        % P converged
        gUpp = gTry;
        PUpp = P; PdotUpp = Pdot; solUpp = sol;
        % Display
        if Display
            fprintf(' gTry = %4.3f \t tConv = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                gTry,P.Time(1),gLow,gUpp);
        end
        
        % Stop Condition
        if (gUpp - gLow - RelTol*gUpp - AbsTol) < AbsTol*1e-3
            break;
        end
    end
end

tTotal = toc(t0);

%% Store Final Result
g = [gLow, gUpp];

if nout > 2
    info = [];
    info.Lower.Gain = gLow;
    info.Lower.P    = PLow;
    info.Lower.Pdot = PdotLow;
    info.Lower.sol  = solLow;
    info.Lower.PowerIterInfo = piterInfo(1:RDEcnt);
    
    info.Upper.Gain = gUpp;
    info.Upper.P    = PUpp;
    info.Upper.Pdot = PdotUpp;
    info.Upper.sol  = solUpp;
    
    info.RDEcnt     = RDEcnt;
    info.TotalTime  = tTotal;
end
end