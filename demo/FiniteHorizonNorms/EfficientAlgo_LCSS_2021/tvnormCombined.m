function [g,d,info] = tvnormCombined(G,varargin)
%% TVNORM1 Computes the time varying system norm using the following combined method:
%
% First, power iterations are performed to compute best lower bound within
% some default tolerance, then we slightly increase the lower bound and
% make sure that RDE has a solution. This is faster than tvnorm which only
% uses bisection appraoch.
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
narginchk(1,4);
nin = nargin;
nout = nargout;

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
[NY,NU] = size(G);

% Obtain default Opt
d  = [];NE = 0;Opt = tvnormOptions;
switch nin       
    case 2
        if isa(varargin{1},'tvnormOptions')
            Opt = varargin{1};
        else
            NE = varargin{1};
        end
    case 3
        NE = varargin{1};
        if isa(varargin{2},'tvnormOptions')
            Opt = varargin{2};
        elseif isa(varargin{2},'tvmat')
            d = varargin{2};
        end
    case 4
        NE = varargin{1};
        d = varargin{2};
        Opt = varargin{3};
end

% DispFlag
DispFlag = isequal(Opt.Display,'on');

% Update cdreOpt based on tvnormoptions
odeOpt = tvodeOptions('OdeOptions',Opt.OdeOptions,'OdeSolver',Opt.OdeSolver);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

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

% Precompute Adjoint System
Ga = tvss(-A',-C(1:NL2,:)',B',D(1:NL2,:)');

%% Initialization
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];

AbsTol     = Opt.AbsTol;
gLow       = Opt.Bounds(1);
gUpp       = Opt.Bounds(2);
GTime      = G.Time;
Nt         = max(20,length(GTime));
x0         = zeros(Nx,1);
RDEcnt     = 0;
MaxPowIter = 500;
piterInfo  = cell(20,1);

if isempty(d)
    d = tvmat(randn(NU,1,Nt),linspace(T0,Tf,Nt));
end
d = d/tvnorm(d);

%% Algorithm
% Start global timing
t0 = tic;

% Main loop
while (gUpp - gLow > AbsTol)
    
    % Power iteration lower bound
    [gPower,d,piterInfo{RDEcnt+1}] = LOCALPowerit(G,Ga,d,CTf,T0,Tf,NL2,x0,MaxPowIter,AbsTol,DispFlag,odeOpt);
    gLow = max(gPower,gLow);
    
    % When you enter the loop, you must try if cleverly chosen bound is
    % an upper bound.
    gUppGuess = gLow+AbsTol;
    gTry = min(gUppGuess,gUpp);
    
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
        if DispFlag
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
        if DispFlag
            fprintf(' gTry = %4.3f \t tConv = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                gTry,P.Time(1),gLow,gUpp);
        end
        
        % Stop Condition
        if (gUpp - gLow - AbsTol) < AbsTol*1e-3
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

function [glb,dwc,info] = LOCALPowerit(G,Ga,d,CTf,T0,Tf,NL2,x0,MaxPowIter,AbsTol,DispFlag,odeOpt)
%% Power iteration lower bound
U      = cell(MaxPowIter,1);
Y      = cell(MaxPowIter,1);
Perf   = zeros(MaxPowIter,1);
U{1}   = d;
t0     = tic;
for i = 1:MaxPowIter
    
    % State Equation
    Y{i} = tvlsim(G,U{i},[T0,Tf],x0,odeOpt);
    
    % Evaluate Performance
    thisY    = Y{i};
    yL2      = thisY(1:NL2,1);
    yE       = tvsubs(thisY(NL2+1:end,1),Tf);
    nyL2     = tvnorm(yL2);
    nyE      = norm(yE);
    Perf(i)  = sqrt(nyE^2 + nyL2^2);
    
    % Display
    if DispFlag
        fprintf(' Iter: %d,\t Perf: %.3f\n',i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxPowIter)
        % Check if performance that we care is stationary
        if Perf(i)-Perf(i-1) < (AbsTol/5)
            break;
        end
    end
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,yL2,[Tf,T0],CTf'*yE,odeOpt);
    
    % Alignment Condition for U
    U{i+1} = U{i+1}/tvnorm(U{i+1});
end

% Final outputs
tTotal = toc(t0);
[glb,mid] = max(Perf(1:i));
dwc = U{mid};
info = struct('TotalTime',tTotal,'Perf',Perf(1:i),'U',U(1:i),'Y',Y(1:i));
end