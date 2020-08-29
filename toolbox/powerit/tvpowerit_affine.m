function [glb,dwc,info] = tvpowerit_affine(G,vbar,NE,d1,x0,Opt)
%% TVPOWERIT_AFFINE
%
% This function performs nominal analysis power iterations for affine
% systems with affine signal provided saperately.
%
% Input argument are as follows:
% G : LTV System
% vbar : fixed additional input
% NE: Number of Euclidean outputs
% d1: starting input for power iterations
% Opt: poweritOptions
%
% The actual input to the system is [d1; vbar]
%
% Outputs are power iteration lower bound, worst-case disturbance and info
% structure.

% Input processing
narginchk(2,6);
nin = nargin;
nout = nargout;
switch nin
    case 2
        NE = 0;d1 = [];x0 = [];Opt = poweritOptions;
    case 3
        d1 = [];x0 = [];Opt = poweritOptions;
    case 4
        x0 = [];Opt = poweritOptions;
    case 5
        Opt = poweritOptions;
end

% Size of the system
[A,B,C,D] = ssdata(G);
NU = size(B,1);
[NY,Nx] = size(C);
NL2 = NY-NE;
Nvbar = size(vbar,1);

% Precompute Adjoint System
Ga = tvss(-A',-C(1:NL2,:)',B',D(1:NL2,:)');

% Horizon and terminal output matrix
[T0,Tf] = getHorizon(G);
CE = tvsubs(C(NL2+1:end,:),Tf);

% Initial condition and input
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(d1)
    d1 = tvmat(randn(NU-Nvbar,1,20),linspace(T0,Tf,20));
else
    Nd = size(d1,1);
    if (NU-Nvbar < Nd)
        error(' Input signal d1 must be of size (total number of plant inputs - total number of affine inputs).');
    end
end

% Options
MaxIter    = Opt.MaxIter;
StopTol    = Opt.StopTol;
OdeSolver  = Opt.OdeSolver;
OdeOptions = Opt.OdeOptions;
Display    = isequal(Opt.Display,'on');
tvopt      = tvodeOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions);

% Memory Allocation
U      = cell(MaxIter,1);
Y      = cell(MaxIter,1);
Perf   = zeros(MaxIter,1);
tComp  = zeros(MaxIter,1);
deltaU = [];

% Initial Input
vtrim = evalt(vbar,d1.Time);
U{1} = [d1/tvnorm(d1); vtrim];

% Start global timing
t0 = tic;

% Power Iterations
for i = 1:MaxIter
    
    % Start timing
    t1 = tic;
    
    % State Equation
    Y{i} = tvlsim(G,U{i},[T0,Tf],x0,tvopt);
    
    % Evaluate Performance
    thisY    = Y{i};
    yL2      = thisY(1:NL2,1);
    yE       = tvsubs(thisY(NL2+1:end,1),Tf);
    nyL2     = tvnorm(yL2);
    nyE      = norm(yE);
    Perf(i)  = sqrt(nyE^2 + nyL2^2);
    
    % Display
    if Display && isempty(deltaU)
        fprintf(dispstr(i,Perf(i),deltaU));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter)
        % Check if performance that we care is stationary
        if abs(Perf(i)-Perf(i-1)) <= StopTol*max(Perf([i i-1]))
            % Check if signals are stationary
            
            [StopTolSatisfied,deltaU] = stopcond(StopTol,U{i},U{i-1});
            if StopTolSatisfied
                tComp(i) = toc(t1);
                break;
            end
        end
    end
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,yL2,[Tf,T0],CE'*yE,tvopt);
    
    % Alignment Condition for U
    nextU = U{i+1};
    vtrim  = evalt(vbar,nextU.Time);
    U{i+1} = [nextU/tvnorm(nextU); vtrim];
    
    % Record single iteration time
    tComp(i) = toc(t1);
end

% Total time
tTotal = toc(t0);

% Final Output
[mPerf,id] = max(Perf(1:i));
glb = mPerf;
dwc = U{id};

if nout > 2
    info = [];
    info.allPerf = Perf(1:i);
    info.allIterComputeTime = tComp(1:i);
    info.TotalTime = tTotal;
    info.TotalIter = i;
end
end