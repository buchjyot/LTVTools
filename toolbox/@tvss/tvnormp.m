function [glb,dwc,info] = tvnormp(G,NE,U1,x0,Opt)
%% TVNORMP
% Runs power iterations to estimate the best lower bound
%
% Input argument are as follows:
% G : LTV System
% NE: Number of Euclidean outputs
% d1: starting input for power iterations
% Opt: poweritOptions
%
% Outputs are power iteration lower bound, worst-case disturbance and info
% structure.

% Input processing
narginchk(1,5);
nin = nargin;
nout = nargout;
switch nin
    case 1
        NE = 0;U1 = [];x0 = [];Opt = poweritOptions;
    case 2
        U1 = [];x0 = [];Opt = poweritOptions;
    case 3
        x0 = [];Opt = poweritOptions;
    case 4
        Opt = poweritOptions;
end

% Precompute Adjoint System
[A,B,C,D] = ssdata(G);
Nu = size(B,2);
[NY,Nx] = size(C);
NL2 = NY-NE;
Ga = tvss(-A',-C(1:NL2,:)',B',D(1:NL2,:)');
[T0,Tf] = getHorizon(G);
CE = tvsubs(C(NL2+1:end,:),Tf);

% Initial condition and input
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(U1)
    U1 = tvmat(randn(Nu,1,20),linspace(T0,Tf,20));
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
deltaU = []; %#ok<NASGU>

% Initial Input
U{1} = U1/tvnorm(U1);

% Display
if Display
    fprintf(' ### Starting Power Iterations:\n');
end

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
    if Display
        fprintf(' Iter: %d,\t Perf: %.3f\n',i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter)
        % Check if performance that we care is stationary
        if Perf(i)-Perf(i-1) <= StopTol*max(Perf([i i-1]))
            % We can check if the input signals are stationary, but the
            % main use of this utility function is to make a "decent" and
            % "quick" guess on the lower bound. Thus, we do not care if the
            % signals are still changing.
            
            % [StopTolSatisfied,deltaU] = stopcond(StopTol,U{i},U{i-1}); if StopTolSatisfied
            tComp(i) = toc(t1);
            break;
            % end
        end
    end
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,yL2,[Tf,T0],CE'*yE,tvopt);
    
    % Alignment Condition for U
    U{i+1} = U{i+1}/tvnorm(U{i+1});
    
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