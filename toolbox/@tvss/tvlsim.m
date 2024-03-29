function [Y,X] = tvlsim(G,U,varargin)
%% TVLSIM Simulate the response of a time-varying system on a specified Tspan
%
% Possible Syntex:
% [Y,X] = tvlsim(G,U,Opt)
% [Y,X] = tvlsim(G,U,T)
% [Y,X] = tvlsim(G,U,T,x0)
% [Y,X] = tvlsim(G,U,T,Opt)
% [Y,X] = tvlsim(G,U,T,x0,Opt)
%
% Where, G is a TVSS
%        [Optional] U is a TVMAT
%        [Optional] T is a time vector that will be supplied to ODE solvers (see
%        documentation of related ODE input arguments for more details)
%        [Optional] x0 is initial condition or terminal condition, default to 0
%        [Optional] Opt is tvodeOptions

%% Input Processing
narginchk(2,5);

% System Data
[A,B,C,D] = ssdata(G);
Nx = order(G);
[Ny,Nu] = size(G);

% Plant horizon
[GT0,GTf] = getHorizon(G);
GTs = G.Ts;
GIM = G.InterpolationMethod;

% Defaults
Tspan = [];x0 = [];Opt = [];

% Input processing
nin = nargin;
switch nin
    case 3
        Vin = varargin{1};
        if isa(Vin,'tvodeOptions')
            Opt = Vin;
        else
            Tspan = Vin;
        end
    case 4
        Tspan = varargin{1};
        Vin = varargin{2};
        if isa(Vin,'tvodeOptions')
            Opt = Vin;
        else
            x0 = Vin;
        end
    case 5
        Tspan = varargin{1};
        x0 = varargin{2};
        Opt = varargin{3};
end

% Defaults
if isempty(Tspan)
    Tspan = [GT0,GTf];
end
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(Opt)
    Opt = tvodeOptions;
end
if isempty(U)
    Tgrid = union(Tspan,G.Time);
    Nt = length(Tgrid);
    if isequal(GTs,0)
        U = tvmat(zeros(Nu,1,Nt),Tgrid,GIM);
    else
        U = tvmat(zeros(Nu,1,Nt),Tgrid,GTs);
    end
end

% Check Sample Time
UTs = U.Ts;
if ~isequal(length(uniquetol([GTs,UTs],sqrt(eps))),1)
    error('Both plant and input must have the same sample time.');
end

%% Simulate System
if isequal(GTs,0)
    [Y,X] = LOCALContinuousTimeSim(G,U,A,B,C,D,Nx,Tspan,x0,Opt);
else
    [Y,X] = LOCALDiscreteTimeSim(A,B,C,D,U,Nx,Ny,GTs,Tspan,x0);
end

function [Y,X] = LOCALDiscreteTimeSim(A,B,C,D,U,Nx,Ny,GTs,Tspan,x0)
%% LOCALDiscreteTimeSim
AData = A.Data;
BData = B.Data;
CData = C.Data;
DData = D.Data;
UData = U.Data;

% Handle Tspan
if length(Tspan) ~= 2
    TimeDiff = uniquetol(diff(Tspan),sqrt(eps));
    if abs(TimeDiff)-GTs > sqrt(eps)
        % When specifying a time vector for time response simulations, the time step must match the sample time of discrete-time models.
        ctrlMsgUtils.error('Control:analysis:rfinputs09');
    end
end

if all(diff(Tspan) > 0)
    % Forward Integration
    SimType = 'Fwd';
    Tgrid = Tspan(1):GTs:Tspan(end);
else
    % Backward Integration
    SimType = 'Bwd';
    Tgrid = Tspan(1):-GTs:Tspan(end);
end

% Time Grid
Nt = length(Tgrid);

% Memory allocation
XData = zeros(Nx,1,Nt);
YData = zeros(Ny,1,Nt);

% Perform integration
switch SimType
    case 'Fwd'
        % Initial Condition
        XData(:,:,1) = x0;
        
        % Forward for Loop
        for i = 1:1:Nt-1
            XData(:,:,i+1) = AData(:,:,i)*XData(:,:,i) + BData(:,:,i)*UData(:,:,i);
            YData(:,:,i)   = CData(:,:,i)*XData(:,:,i) + DData(:,:,i)*UData(:,:,i);
        end
        YData(:,:,Nt) = CData(:,:,Nt)*XData(:,:,Nt) + DData(:,:,Nt)*UData(:,:,Nt);
        
        % Final Output
        X = tvmat(XData,Tgrid,GTs);
        Y = tvmat(YData,Tgrid,GTs);
        
    case 'Bwd'
        % Final Condition
        XData(:,:,Nt) = x0;
        
        % Backward for Loop
        for i = Nt:-1:2
            XData(:,:,i-1) = AData(:,:,i-1)*XData(:,:,i) + BData(:,:,i-1)*UData(:,:,i-1);
            YData(:,:,i)   = CData(:,:,i)*XData(:,:,i) + DData(:,:,i)*UData(:,:,i);
        end
        YData(:,:,1)   = CData(:,:,1)*XData(:,:,i) + DData(:,:,i)*UData(:,:,i);
        
        % Final Output
        X = tvmat(XData,flip(Tgrid),GTs);
        Y = tvmat(YData,flip(Tgrid),GTs);
end

function [Y,X] = LOCALContinuousTimeSim(G,U,A,B,C,D,Nx,Tspan,x0,Opt)
%% LOCALContinuousTimeSim
% If Nx = 0 then return Y = D*U;
if isequal(Nx,0)
    Y = evalt(D,U.Time)*U;
    X = [];
    return;
end

% Read Varargins
OdeSolverFh = str2func(Opt.OdeSolver);
OdeOpt = Opt.OdeOptions;

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
[t,x] = OdeSolverFh(odefh,Tspan,x0,OdeOpt);

%% Construct Outputs

% State vector X
% TVMAT must have increasing time grid points
if all(diff(Tspan)>0)
    X = tvmat(x,t);
else
    X = tvmat(flip(x),flip(t));
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