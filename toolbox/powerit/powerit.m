function varargout = powerit(mySystem,varargin)
%% POWERIT
% This function finds the worst-case input that maximizes the induced gain
%
%   OBJ = POWERIT(SYS)  Finds the input that maximizes the induced L2-to-L2
%   gain of the system by simulating the system on default "fixed" horizon.
%
%       SYS could be one of the following:
%       1) LTI state-space systems (SS object)
%       2) time-varying state-space (TVSS object)
%       3) polynomial dynamic system (POLYSYS object)
%       4) simulink model (character vector specifying name of the model)
%       5) m-file S-function (specified using function_handle)
%
%       The related default fixed horizon is obtained using the following rule:
%       1) default "fixed" horizon [T0,Tf] = [0,1] sec
%       2) [T0,Tf] = getHorizon(SYS);
%       3) default "fixed" horizon [T0,Tf] = [0,1] sec
%       4) simulate the model on pre-configured horizon and determine [T0,Tf]
%       5) default "fixed" horizon [T0,Tf] = [0,1] sec
%
%   OBJ = POWERIT(SYS,tUser) User can provide the time vector/horizon
%   "tUser" on which analysis can be performed. tUser must be numeric and
%   of length > 1. It must be in the form of tUser = [T0,TF] or tUser =
%   [T0,T1,...,TF]. This option overrides the default horizon.
%
%   [OBJ,dWC,INFO] = POWERIT(...) In addition to objective, user can get
%   worst-case disturbance signal dWC. Moreover, INFO provides additional
%   info such as TotalTime to run all the iterations and worst-case
%   states/outputs.
%
%   [...] = POWERIT(SYS,tUser,U1,x0,Opt) In addition to other arguments,
%   user can provide initial input U1 as TVMAT, x0 as a vector of (Nx by 1)
%   initial conditions and Opt specified through poweritOptions.
%
%       NOTE: These additional inputs (U1,x0,Opt) override default
%       configuration. x0 can not be specified for simulink models, as we
%       assume that initial conditions must be part of the simulink model.
%
%   [...] = POWERIT(SYS,tUser,U1,x0,SIGSPEC,Opt) User can provide a
%   poweritSignalSpec
%
%       Default values are (NE,Np,Nw,Nv) = (0,0,0,0) and UncParamRange =
%       [0,0]. If system has Ny outputs, and user specifies NE = Ny, then
%       (NE,Np,Nw,Nv) = (Ny,0,0,0) computes nominal induced L2-to-Euclidean
%       gain of the system. In other cases (i.e. NE < Ny), we compute the
%       gain from L2 inputs to a mixture of L2 and Euclidean penalty on
%       outputs.
%
%   See also poweritOptions.

%% Input Processing
narginchk(1,8);
nin  = nargin;
nout = nargout;

%% System Type
systemClass = class(mySystem);
switch systemClass
    case {'ss','tvss','polysys'}
        systemType = systemClass;
        
    case 'char'
        systemType = 'model';
        
    case 'function_handle'
        systemType = 'sfunc';
        
    otherwise
        error('Invalid input system.');
end

%% User Time Grid
% First input must be the time grid
tUser = [];
if nin > 1
    tUser = varargin{1};
    varargin = varargin(2:end);
else
    switch systemType
        case {'ss','polysys','sfunc'}
            tUser = [0,1];
            
        case 'tvss'
            [T0,Tf] = getHorizon(mySystem);
            tUser = [T0,Tf];
            
        case 'model'
            % This assumes that simulink model will stop at some fixed time.
            tnom = sim(mySystem);
            tUser = [tnom(1),tnom(end)];
    end
end

% isValid
isValidTimeVec = isnumeric(tUser)  &&  (length(tUser) > 1);
if ~isValidTimeVec
    error('Time vector must be an array of the form T = [T0,TF] or T = [T0,T1,...,TF].');
end

% Lift up the "ss" object to TVMAT
if isequal(systemType,'ss')
    mySystem = tvss(mySystem,tUser);
    systemType = 'tvss';
end

%% Options
id = cellfun(@(x) isa(x,'poweritOptions'),varargin);
V1 = varargin(id);
if isempty(V1)
    pOpt = poweritOptions;
else
    pOpt = V1{1};
end
varargin = varargin(~id);

%% IO Sizes
% Get System Sizes
[Nx,NOut,NIn] = wcgetsizes(mySystem,systemType);
id = cellfun(@(x) isa(x,'poweritSignalSpec'),varargin);
V2 = varargin(id);
if isempty(V2)
    pSpec = poweritSignalSpec;
else
    pSpec = V2{1};
    varargin = varargin(~id);
end
warn_state = warning('off','MATLAB:structOnObject');
pSpec = struct(pSpec);
warning(warn_state);
fnames = fieldnames(pSpec);
Nfield = length(fnames);
for k = 1:Nfield
    eval([fnames{k} '= pSpec.(fnames{k});']);
end

Nu = NIn-Np-Nw;
pSpec.Nx = Nx;
pSpec.NOut = NOut;
pSpec.NIn = NIn;
pSpec.Nu = Nu;
pSpec.NL2 = NOut-NE-Nv;

%% Initial Conditions
% Validate InitialCondCostMat
if isa(InitialCondCostMat,'function_handle')
    E0 = InitialCondCostMat(Nx);
else
    % double Nx by Nx
    E0 = InitialCondCostMat;
    nr = size(E0,1);
    if ~isequal(nr,Nx)
        error('InitialCondCostMat must be square matrix of size Nx by Nx.')
    end
end

% Read initial conditions
fixedIC = isequal(InitialConditions,'fixed');
id = cellfun(@(x) isa(x,'double'),varargin);
V3 = varargin(id);
zeroIC = [];
if isempty(V3)
    if fixedIC
        % Specify default initial conditions if x0 is empty
        switch systemType
            case {'tvss','polysys'}
                x0 = zeros(Nx,1);
            case 'sfunc'
                [~,x0] = feval(mySystem,[],[],[],0);
            case 'model'
                x0 = [];
        end
    else
        % Free initial condition for states
        zeroIC = isinf(diag(E0));
        E0 = E0(~zeroIC,~zeroIC);
        x0 = randn(size(E0,1),1);
        x0 = x0/(x0'*E0*x0);
    end
else
    x0 = V3{1};
    % Make sure the length of initial condition is appropriate
    if ~isequal(length(x0(:)),Nx)
        error('Initial conditions vector must be of length Nx.');
    end
end
varargin = varargin(~id);

% Update InitialCondCostMat
pSpec.zeroIC = zeroIC;
pSpec.InitialCondCostMat = E0;

% If user specified initial conditions then error out, because in
% this case user can not explicitly specify initial conditions, it
% has to be part of the simulink model.
if ~isempty(x0) && isequal(systemType,'model')
    error('For Simulink model initial conditions must be fixed and part of the model.');
end

%% Initial Input
id = cellfun(@(x) isa(x,'tvmat'),varargin);
V4 = varargin(id);
if isempty(V4)
    InitialInput = pSpec.InitialInput;
    IntialParam  = sum(pSpec.UncParamRange,2)/2;
    fcn = str2func(InitialInput);
    Nt  = length(tUser);
    if isequal(Nt,2) && (isequal(InitialInput,'randn') || isequal(InitialInput,'rand'))
        % In this case user just provided the horizon for analysis. We
        % should define the input with some more samples to have more
        % exploration. We choose sample time as 5% of the horizon.
        Td = tUser(end)*(5/100);
        U1Time = unique([tUser(1) tUser(1):Td:tUser(end) tUser(end)]);
        Nt = length(U1Time);
        U1 = tvmat([fcn(Nw+Nu,1,Nt); IntialParam.*ones(Np,1,Nt)],U1Time);
    else
        % In this case, user provided a grid, so just use it to define the
        % initial input.
        U1  = tvmat([fcn(Nw+Nu,1,Nt); IntialParam.*ones(Np,1,Nt)],tUser);
    end
else
    U1 = V4{1};
    [Nr,Nc] = size(U1);
    if ~isequal(Nr,Nw+Nu+Np) || ~isequal(Nc,1)
        error('Initial input dimentions must agree with input dimentions of the model.');
    end
end
U1 = U1*(InputL2Norm/tvnorm(U1));

%% Call Solver
% Suffix ic denotes that these solvers also consider constraint on intial
% conditions.
isUncertain = (Np||Nv||Nw);
if isUncertain
    % Robust Analysis
    solverfh = @uwcsolver_default;
else
    % Nominal Analysis
    if fixedIC
        solver_case = 1;
        switch solver_case
            case 1
                solverfh = @wcsolver_default;
            case 2
                solverfh = @wcsolver_default_normalize_twice;
            case 3
                solverfh = @wcsolver_yFminusy0_terminalconstraint;
        end
    else
        % solverfh = @wcsolver_combined_initialcondition_inputnorm_constraint;
        solverfh = @wcsolver_separate_initialcondition_inputnorm_constraint;
    end
end
[varargout{1:max(nout,1)}] = feval(solverfh,mySystem,tUser,U1,x0(:),pSpec,pOpt,systemType);
end

function [nStates,nOutputs,nInputs] = wcgetsizes(mySystem,systemType)
%% WCGETSIZES
switch systemType
    case 'tvss'
        nStates = order(mySystem);
        [nOutputs,nInputs] = size(mySystem);
        
    case 'model'
        sizeInfo = feval(mySystem,[],[],[],'sizes');
        nStates  = sizeInfo(1);  % Number of continuous states
        nOutputs = sizeInfo(3);  % Number of outputs
        nInputs  = sizeInfo(4);  % Number of inputs
        % XXX - simulink models with discrete states are not supported.
        
    case 'sfunc'
        sizeInfo = feval(mySystem,[],[],[],0);
        nStates  = sizeInfo(1);  % Number of continuous states
        nOutputs = sizeInfo(3);  % Number of outputs
        nInputs  = sizeInfo(4);  % Number of inputs
        
    case 'polysys'
        nStates = length( mySystem.states );
        nInputs = length( mySystem.inputs );
        nOutputs = size( mySystem.orMap, 1 );
end
end