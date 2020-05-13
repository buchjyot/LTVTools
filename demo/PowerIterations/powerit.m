function [glb,dwc,info] = powerit(G,varargin)
%% POWERIT
%
% Input:
% G     = must be TVSS or Simulink model on a finite horizon
% Opt   = options for power iterations
%
% Output:
% glb   = lower bound on finite horizon performance
% dwc   = final worst-case disturbance
% info  = additional info
%
% This is a root level function that makes a call to utility function

%% Input Processing
narginchk(1,2);

% Read Options
Opt = varargin{cellfun(@(x) isa(x,'poweritOptions'),varargin)};
if isempty(Opt), Opt = poweritOptions; end

% Options
MaxIter         = Opt.MaxIter;
StopTol         = Opt.StopTol;
OdeSolver       = Opt.OdeSolver;
OdeOptions      = Opt.OdeOptions;
StepSize        = Opt.StepSize;
Objective       = Opt.Objective;

% Options
Display         = isequal(Opt.Display,'on');
FixedStepSize   = isa(StepSize,'double');
tvsopt1         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','ForwardInTime','StepSize',StepSize);
tvsopt2         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','BackwardInTime','StepSize',StepSize);

%% Fixed Horizon and I/O Sizes
switch class(G)
    case 'char'
        % We assume that the simulink model is configured such that it will
        % stop at some finite final time Tf.
        load_system(G);
        [tnom,xnom,ynom] = sim(G); %#ok<*ASGLU>
        T0 = tnom(1);
        Tf = tnom(end);
        
        % I/O Sizes
        io  = getlinio(G);
        Nu = evalin('base','size(simin,2)-1');
        
    case 'tvss'
        [~,Nu]  = size(G);
        [T0,Tf] = getHorizon(G);
end

% StepSize
if FixedStepSize
    Tgrid = T0:StepSize:Tf;
else
    Tgrid = [T0,Tf];
end

%% Memory Allocation
U    = cell(MaxIter,1);
Y    = cell(MaxIter,1);
Perf = zeros(MaxIter,1);

%% Define Functions
% Predefine functions to avoid overhead of if/else statement in the loop
switch class(G)
    case 'char'
        % Simulink Model
        SimulateStateEquations      = @(Ui) LOCALSim(G,Tgrid,Ui);
        SimulateCostateEquations    = @(Yi) tvlsim(LOCALAdjoint(G,io,Yi.Time),Yi,tvsopt2);
        U{1}                        = validateInitialInput(Opt,tnom,Nu);
        
    case 'tvss'
        % LTV Model
        SimulateStateEquations      = @(Ui) tvlsim(G,Ui,tvsopt1);
        Ga                          = G';
        SimulateCostateEquations    = @(Yi) tvlsim(Ga,Yi,tvsopt2);
        U{1}                        = validateInitialInput(Opt,G.Time,Nu);
end

% Display String & Performance
DisplayProgress = @(id,pf) fprintf('Iter: %d, Performance: %.3f\n',id,pf);
switch Objective
    case 'L2toL2'
        EvaluatePerformance = @(in) tvnorm(in);
    case 'L2toE'
        EvaluatePerformance = @(in) norm(tvsubs(in,in.Time(end)));
end

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    Y{i} = SimulateStateEquations(U{i});
    
    % Alignment Condition for Y
    Perf(i) = EvaluatePerformance(Y{i});
    Y{i} = Y{i}/Perf(i);
    
    % Display
    if Display
        DisplayProgress(i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter-1)
        % Performance that we care is stationary
        if abs(Perf(i)-Perf(i-1)) <= StopTol*max(Perf([i i-1]))
            % Performance may reach stationary but input signal may not
            if FixedStepSize
                if tvnorm(U{i}-U{i-1}) <= StopTol
                    break;
                end
            else
                % For variable stepsize, we will have two different time
                % grids, thus for subtraction we need consistent time grid.
                tUnion = union(U{i}.Time,U{i-1}.Time);
                [Ui,Uim1] = evalt(U{i},U{i-1},tUnion);
                if tvnorm(Ui-Uim1) <= StopTol
                    break;
                end
            end
        end
    end
    
    % Costate Equation
    U{i+1} = SimulateCostateEquations(Y{i});
    
    % Alignment Condition for U
    U{i+1} = U{i+1}/tvnorm(U{i+1});
end
tcomp = toc(t1);

%% Final Output
glb = Perf(i);
dwc = U{i};

info = [];
info.TotalIter  = i;
info.TotalTime  = tcomp;
end

function Ga = LOCALAdjoint(model,io,tSnapShot)
%% Linearize Dynamics and Construct Adjoint System
[linsys,linop]          = linearize(model,io,tSnapShot);
Glin                    = tvss(linsys,tSnapShot);
[dfdx,dfdu,dgdx,dgdu]   = ssdata(Glin);
Ga                      = tvss(-dfdx',-dgdx',dfdu',dgdu');
end

function Yi = LOCALSim(model,Tgrid,Ui)
%% Simulate Simulink Model
assignin('base','simin',[Ui.Time reshapedata(Ui)]);
[t,x,y] = sim(model,Tgrid);
Yi = tvmat(y,t);
end