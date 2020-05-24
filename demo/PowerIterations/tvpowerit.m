function [glb,dwc,info] = tvpowerit(G,varargin)
%% TVPOWERIT
%
% Input:
% G     = must be TVSS on a finite horizon
% NE    = [Optional] number of Euclidean outputs
% x0    = [Optional] initial conditions
% U1    = [Optional] initial input for power iterations
% Opt   = [Optional] options for power iterations (tvpoweritOptions)
%
% Output:
% glb   = lower bound on finite horizon performance
% dwc   = final worst-case disturbance
% info  = additional info
%
% Possible Syntex:
% [glb,dwc,info] = tvpowerit(G);
% [glb,dwc,info] = tvpowerit(G,NE);
% [glb,dwc,info] = tvpowerit(G,NE,x0);
% [glb,dwc,info] = tvpowerit(G,NE,x0,U1);
% [glb,dwc,info] = tvpowerit(G,NE,x0,U1,Opt);

%% Input Processing
narginchk(1,5);

% Input Parser
[NE,x0,U1,Opt] = LOCALInputParser(varargin);

% Options
MaxIter          = Opt.MaxIter;
StopTol          = Opt.StopTol;
OdeSolver        = Opt.OdeSolver;
OdeOptions       = Opt.OdeOptions;
StepSize         = Opt.StepSize;
Display          = isequal(Opt.Display,'on');
FixedStepSize    = isa(StepSize,'double');
tvsopt1          = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','ForwardInTime','StepSize',StepSize);
tvsopt2          = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','BackwardInTime','StepSize',StepSize);
StopTolSatisfied = false;

%% System Properties
Nx      = order(G);
[Ny,Nu] = size(G);
[~,Tf]  = getHorizon(G);
NL2     = Ny-NE;

%% Memory Allocation
U    = cell(MaxIter,1);
Y    = cell(MaxIter,1);
X    = cell(MaxIter,1);
W    = cell(MaxIter,1);
Perf = zeros(MaxIter,1);

%% Define Functions
% Display String & Performance
if NE == 0
    % L2toL2
    EvaluatePerformance = @(Yi,Xi) tvnorm(Yi);
    lamTf               = @(Xi) zeros(Nx,1);
    getCostateInputs    = @(Yi,pf) Yi/pf;
elseif NE == Ny
    % L2toE
    F                   = LOCALGetTerminalCostMatrix(G,NL2,Ny,Tf);
    EvaluatePerformance = @(Yi,Xi) sqrt(tvsubs(Xi,Tf)'*F*tvsubs(Xi,Tf));
    lamTf               = @(Xi) F*tvsubs(Xi,Tf);
    getCostateInputs    = @(Yi,pf) evalt(tvmat(zeros(NE,1)),Yi.Time);
else
    % Mixed (L2 and Euclidean) Cost
    F                   = LOCALGetTerminalCostMatrix(G,NL2,Ny,Tf);
    EvaluatePerformance = @(Yi,Xi) sqrt(tvsubs(Xi,Tf)'*F*tvsubs(Xi,Tf) + tvnorm(Yi(1:NL2))^2);
    lamTf               = @(Xi) F*tvsubs(Xi,Tf);
    getCostateInputs    = @(Yi,pf) [Yi(1:NL2)/tvnorm(Yi(1:NL2)); zeros(NE,1)];
end

% Adjoint System
Ga = G';

% Stop Condition
if FixedStepSize
    StopCondition = @(Ui,Uim1) tvnorm(Ui-Uim1) <= StopTol;
else
    % If stepsize is not fixed then we need to evaluate Ui and Uiminus1 on
    % the same time grid in order to compare them for stopping condition.
    UnionTimeGrid = @(Ui,Uim1) union(Ui.Time,Uim1.Time);
    StopCondition = @(Ui,Uim1) tvnorm(evalt(Ui,UnionTimeGrid(Ui,Uim1))-evalt(Uim1,UnionTimeGrid(Ui,Uim1))) <= StopTol;
end

%% Validate Initial Input
if isempty(U1)
    % Default input is "randn"
    Nt  = length(G.Time);
    U1 = tvmat(randn(Nu,1,Nt),G.Time);
else
    [Nr,Nc] = size(U1);
    if ~isequal(Nr,Nu) || ~isequal(Nc,1)
        error('Initial input dimentions must agree with input dimentions of the model.');
    end
end

% Normalize Initial Input
U{1} = U1/tvnorm(U1);

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    [Y{i},X{i}] = tvlsim(G,U{i},x0,tvsopt1);
    
    % Evaluate Performance
    Perf(i) = EvaluatePerformance(Y{i},X{i});
    
    % Display
    if Display
        fprintf(' Iter: %d, Performance: %.3f\n',i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter-1)
        % Check if performance that we care is stationary
        if abs(Perf(i)-Perf(i-1)) <= StopTol*max(Perf([i i-1]))
            % Check if input signal is stationary
            if StopCondition(U{i},U{i-1})
                StopTolSatisfied = true;
                break;
            end
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
tTotal = toc(t1);

%% Final Output
[mPerf,id] = max(Perf(1:i));
glb = mPerf;
dwc = U{id};
if Display
    if isequal(i,MaxIter) && ~StopTolSatisfied
        fprintf(' Maximum number of iterations reached.\n');
    end
    if StopTolSatisfied
        fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
    end
end

info = [];
info.TotalIter = i;
info.TotalTime = tTotal;
info.allPerf = Perf(1:i);
end

function [NE,x0,U1,Opt] = LOCALInputParser(V)
%% LOCALInputParser

% Defaults
NE  = 0;
x0  = [];
U1  = [];
Opt = tvpoweritOptions;

% Identify distinct arguments
V1 = V(cellfun(@(x) isa(x,'tvmat'),V));
if ~isempty(V1)
    U1 = V1{1};
end

V2 = V(cellfun(@(x) isa(x,'tvpoweritOptions'),V));
if ~isempty(V2)
    Opt = V2{1};
end

% x0 and NE are both double, first argument has to be NE and next double
% argument must be x0.
V3d  = V(cellfun(@(x) isa(x,'double'),V));
if isempty(V3d)
    return;
elseif isequal(length(V3d),2)
    NE = V3d{1};
    x0 = V3d{2};
else
    NE = V3d{1};
end
end

function F = LOCALGetTerminalCostMatrix(G,NL2,Ny,Tf)
%% LOCALGetTerminalCostMatrix
% Extract Euclidean part
[~,~,CE,DE] = ssdata(tvss(G.Data(NL2+1:Ny,:),G.Time,G.InterpolationMethod));
[CTf,DTf] = tvsubs(CE,DE,Tf);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm');
end

% Terminal Cost Matrix
F = CTf'*CTf;
end