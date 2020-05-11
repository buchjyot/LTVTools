function [Y,X] = tvlsim(G,U,varargin)
%% TVLSIM Simulate the response of a time-varying system
%
% [Y,X] = tvlsim(G,U)
% [Y,X] = tvlsim(G,U,x0)
% [Y,X] = tvlsim(G,U,Opt)
% [Y,X] = tvlsim(G,U,x0,Opt)
%
% Where, G is a TVSS
%        U is a TVMAT
%        [Optional] x0 is initial condition or terminal condition, default to 0
%        [Optional] Opt is tvlsimOptions

%% Input Processing
narginchk(2,4);
nin = nargin;

% System Data
[A,B,C,D] = ssdata(G);
Nx = order(G);

% If Nx = 0 then return Y = D*U;
if isequal(Nx,0)
    Y = evalt(D,U.Time)*U;
    X = [];
    return;
end

% Read Options
Opt = [];x0 = [];
switch nin
    case 3
        if isa(varargin{1},'tvlsimOptions')
            Opt = varargin{1};
        else
            x0 = varargin{1};
        end
    case 4
        x0 = varargin{1};
        Opt = varargin{2};
end

% Assume zero boundary conditions if user did not provide any
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(Opt)
    Opt = tvlsimOptions;
end
OdeSolverFh = str2func(Opt.OdeSolver);
OdeOpt = Opt.OdeOptions;

% Make sure the system (G) is defined on the given input (U) horizon
ltvutil.verifyFH(G,U);
[GT0,GTf] = getHorizon(G);
[UT0,UTf] = getHorizon(U);
% NOTE: The following code allows the simulation to happen even if the
% input is only defined on a part of the horizon.
if any(UT0>GTf || UTf>GTf || UT0<GT0 || UTf<GT0)
    error('Input U must be defined on the subset of the system horizon.');
end

%% Tgrid
% Integration time grid
StepSize = Opt.StepSize;
switch StepSize
    case 'Default'
        % Similar to MATLAB command lsim, this 'Default' step size assumes
        % integration to be performed on the time grid specified by U.Time
        % The returned output is also in the same time grid as U.Time.
        Tgrid = U.Time;
        
    case 'Auto'
        % Time grid is automatically determined by ODE solver, we only
        % provide the span of integration
        Tgrid = [GT0,GTf];
        
    otherwise
        % This means StepSize is "Fixed" to the specified double value
        Tgrid = GT0:StepSize:GTf;
end

% MATLAB command lsim does not allow backward simulation, whereas tvlsim
% allows backward simulation if simulation "Type" option is specified as "BackwardInTime".
if isequal(Opt.Type,'BackwardInTime')
    Tgrid = flip(Tgrid);
end

%% Simulate System
% Specify ODE
UseLinearInterp = isequal(G.InterpolationMethod,'Linear') && isequal(U.InterpolationMethod,'Linear');
if UseLinearInterp
    AB = [A B];
    ABData = AB.Data;
    ABTime = AB.Time;
    ABDiff = diff(ABTime);
    NABTime = numel(ABTime);
    
    UData = U.Data;
    UTime = U.Time;
    UDiff = diff(UTime);
    NUTime = numel(UTime);
    
    odefh = @(t,x) LOCALderivLinInterp(t,x,ABTime,NABTime,ABDiff,ABData,...
        UTime,NUTime,UDiff,UData);
else
    odefh = @(t,x) LOCALderiv(t,x,A,B,U);
end

% Solve ODE
[t,x] = OdeSolverFh(odefh,Tgrid,x0,OdeOpt);

%% Construct Outputs

% State vector X
% TVMAT must have increasing time grid points
if isequal(Opt.Type,'BackwardInTime')
    X = tvmat(flip(x),flip(t));
else
    X = tvmat(x,t);
end

[U1,C,D] = evalt(U,C,D,X.Time);
Y = C*X+D*U1;


function xdot = LOCALderiv(t,x,A,B,U)
%% LOCALderiv
[At,Bt,ut] = tvsubs(A,B,U,t);
xdot = At*x + Bt*ut;


function xdot = LOCALderivLinInterp(t,x,ABTime,NABTime,ABDiff,ABData,...
    UTime,NUTime,UDiff,UData)
%% LOCALderivLinInterp
% Eval AB(t)
if t<=ABTime(end) && t>=ABTime(1)
    % Within the horizon inclusive of boundary points
    [k,alpha] = LOCALfindslotalpha(NABTime,ABTime,t,ABDiff);
    if alpha==0
        ABt = ABData(:,:,k);
    else
        ABt = (1-alpha)*ABData(:,:,k) + alpha*ABData(:,:,k+1);
    end
else
    % Extrapolate to zeros if outside horizon
    [nr,nc] = size(ABData(:,:,1));
    ABt = zeros(nr,nc);
end

% Eval U(t)
if t<=UTime(end) && t>=UTime(1)
    % Within the horizon inclusive of boundary points
    [k,alpha] = LOCALfindslotalpha(NUTime,UTime,t,UDiff);
    if alpha==0
        Ut = UData(:,:,k);
    else
        Ut = (1-alpha)*UData(:,:,k) + alpha*UData(:,:,k+1);
    end
else
    % Extrapolate to zeros if outside horizon
    [nr,nc] = size(UData(:,:,1));
    Ut = zeros(nr,nc);
end

xdot = ABt*[x; Ut];


function [k,alpha] = LOCALfindslotalpha(N,vec,val,dvec)
%% LOCALfindslotalpha
% N integer
% vec 1-by-N (or N-by-1), sorted
% val 1-by-1
% dvec = diff(vec)

k = max(find(val>=vec));   %#ok<MXFND> % don't follow advice - it is slower.
if ~isempty(k)
    if k<N
        alpha = (val - vec(k))/dvec(k);
    else
        alpha = 0;
    end
else
    k = 1;
    alpha = 0;
end
