function [glb,dwc,info] = slpowerit(model,varargin)
%% SLPOWERIT
%
% Input:
% model = must be a Simulink model set up to stop at finite horizon
% U1    = [Optional] initial input
% Opt   = [Optional] options for power iterations [slpoweritOptions]
%
% Output:
% glb   = lower bound on finite horizon performance
% dwc   = final worst-case disturbance
% info  = additional info
%
% This function assumes that final time is fixed.

%% Input Processing
narginchk(1,2);

% Model has to be simulink model
if ~ischar(model)
    error('Input model must be character vector specifying Simulink model.');
end

% Input Parser
[U1,Opt] = LOCALInputParser(varargin);

% Options
MaxIter          = Opt.MaxIter;
StopTol          = Opt.StopTol;
OdeSolver        = Opt.OdeSolver;
OdeOptions       = Opt.OdeOptions;
StepSize         = Opt.StepSize;
Objective        = Opt.Objective;
LinOpt           = Opt.LinOpt;
InputL2Norm      = Opt.InputL2Norm;
Display          = isequal(Opt.Display,'on');
FixedStepSize    = isa(StepSize,'double');
tvsopt2          = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','BackwardInTime','StepSize',StepSize);
StopTolSatisfied = false;

%% Fixed Horizon and I/O Sizes
% We assume that the simulink model is configured such that it will
% stop at some "fixed" finite final time Tf.
load_system(model);
[tnom,xnom,ynom] = sim(model); %#ok<*ASGLU>

% I/O Sizes
io  = getlinio(model);
Nu = evalin('base','size(simin,2)-1');
Nx = size(xnom,2);
Ny = size(ynom,2);

% Nominal Horizon
Tf0 = tnom(end);
if Display
    fprintf(' Nominal Horizon: %.3f\n',Tf0);
end

%% Memory Allocation
U    = cell(MaxIter,1);
Y    = cell(MaxIter,1);
X    = cell(MaxIter,1);
W    = cell(MaxIter,1);
Perf = zeros(MaxIter,1);

%% Validate Initial Input
if isempty(U1)
    % Default input is "randn"
    Nt = length(tnom);
    U1 = tvmat(randn(Nu,1,Nt),tnom);
else
    [Nr,Nc] = size(U1);
    if ~isequal(Nr,Nu) || ~isequal(Nc,1)
        error('Initial input dimentions must agree with input dimentions of the model.');
    end
end

% Normalize Initial Input
U{1} = U1*(InputL2Norm/tvnorm(U1));

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    assignin('base','simin',[U{i}.Time reshapedata(U{i})]);
    if FixedStepSize
        [Tgrid1,x,y] = sim(model,'FixedStep',num2str(StepSize));
        % [Tgrid1,x,y] = sim(model,Tf0,'FixedStep',num2str(StepSize));
    else
        [Tgrid1,x,y] = sim(model);
        % [Tgrid1,x,y] = sim(model,Tf0);
    end
    Y{i} = tvmat(y,Tgrid1);
    X{i} = tvmat(x,Tgrid1);
    
    % Compute Costate Dynamics using Snapshot linearization
    tSnapShot             = Y{i}.Time;
    [linsys,linop]        = linearize(model,io,tSnapShot,LinOpt);
    Tgrid2                = linsys.SamplingGrid.Time;
    
    % If the model is configured to stop at some terminal condition then
    % the linsys may have few less time grid points (at the end) than
    % requested by tSnapShot. The following code assumes that both time
    % grids are identical.
    Glin                  = tvss(linsys,Tgrid2);
    [dfdx,dfdu,dgdx,dgdu] = ssdata(Glin);
    
    % Construct Adjoint System
    Ga = tvss(-dfdx',-dgdx',dfdu',dgdu');
    
    % If the model is configured to stop at some terminal condition then
    % the linsys may have few less time grid points (at the end) than
    % requested. If the final times are different then choose the
    % linearized system time grid instead.
    if ~isequal(Tgrid1(end),Tgrid2(end))
        Y{i} = evalt(Y{i},Tgrid2);
        X{i} = evalt(X{i},Tgrid2);
    end
    
    % Alignment Condition for Y
    Tf1 = Tgrid2(end);
    switch Objective
        case 'L2toL2'
            Perf(i) = tvnorm(Y{i})/InputL2Norm;
            lamTf   = zeros(Nx,1);
            W{i}    = Y{i}/Perf(i);
        case 'L2toE'
            Perf(i) = sqrt(tvsubs(Y{i},Tf1)'*tvsubs(Y{i},Tf1))/InputL2Norm;
            lamTf   = tvsubs(dgdx,Tf1)'*tvsubs(dgdx,Tf1)*tvsubs(X{i},Tf1);
            W{i}    = evalt(tvmat(zeros(Ny,1)),Tgrid2);
    end
    
    % Display
    if Display
        fprintf(' Iter: %d, Horizon: %.3f, Performance: %.3f\n',i,Tf1,Perf(i));
    end
    
    % Stopping Condition
    if (i > 2) && (i < MaxIter-1)
        % Performance that we care is stationary
        if abs(Perf(i)-Perf(i-1)) <= StopTol*max(Perf([i i-1]))
            % Performance may reach stationary but input signal may not
            if StopCondition(U{i},U{i-1},FixedStepSize,StopTol)
                StopTolSatisfied = true;
                break;
            end
        end
    end
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,W{i},lamTf,tvsopt2);
    
    % Alignment Condition for U
    U{i+1} = U{i+1}*(InputL2Norm/tvnorm(U{i+1}));
end
tcomp = toc(t1);

%% Final Output
[mPerf,id] = max(Perf(2:i));
glb = mPerf;
dwc = U{id+1};
if Display
    if isequal(i,MaxIter) && ~StopTolSatisfied
        fprintf(' Maximum number of iterations reached.\n');
    end
    if StopTolSatisfied
        fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
    end
end

info           = [];
info.TotalIter = i;
info.TotalTime = tcomp;
info.allPerf   = Perf(1:i);
end

function out = StopCondition(Ui,Uim1,FixedStepSize,StopTol)
%% Stop Condition
if FixedStepSize
    out = tvnorm(Ui-Uim1) <= StopTol;
else
    tUnion = union(Ui.Time,Uim1.Time);
    Tfmin = min(Ui.Time(end),Uim1.Time(end));
    tUnion(tUnion > Tfmin) = [];
    out = tvnorm(evalt(Ui,tUnion)-evalt(Uim1,tUnion)) <= StopTol;
end
end

function [U1,Opt] = LOCALInputParser(V)
%% LOCALInputParser

% Defaults
U1  = [];
Opt = slpoweritOptions;

% Identify distinct arguments
V1 = V(cellfun(@(x) isa(x,'tvmat'),V));
if ~isempty(V1)
    U1 = V1{1};
end

V2 = V(cellfun(@(x) isa(x,'slpoweritOptions'),V));
if ~isempty(V2)
    Opt = V2{1};
end
end