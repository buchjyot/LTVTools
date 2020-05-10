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

Nx = order(G);
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

% If empty then set to default
if isempty(x0)
    % Assume zero boundary conditions if user did not provided any i.e.
    % tvlsim(G,U,[],Opt)
    x0 = zeros(Nx,1);
end
if isempty(Opt)
    Opt = tvlsimOptions;
end
OdeSolverFh = str2func(Opt.OdeSolver);
OdeOpt = Opt.OdeOptions;

% Make sure the system (G) is defined on the given input (U) horizon 
[A,B,C,D] = ssdata(G);
ltvutil.verifyFH(G,U);
[GT0,GTf] = getHorizon(G);
[UT0,UTf] = getHorizon(U);
% NOTE: The following code allows the simulation to happen in the subset of
% the plant horizon, even if the input is only defined on part of the
% horizon. Where the input is not defined, the "tvsubs" function used in
% the LOCALderiv will extrapolate to zero.
if any(UT0>GTf || UTf>GTf || UT0<GT0 || UTf<GT0)
    error('Input U must be defined on the subset of the system horizon.');
end

% Tspan
if isequal(Opt.Type,'ForwardInTime')
    Tspan = [GT0 GTf];
else
    Tspan = [GTf GT0];
end

%% Simulate System
% Solve ODE
odefh = @(t,x) LOCALderiv(t,x,A,B,U);
[t,x] = OdeSolverFh(odefh,Tspan,x0,OdeOpt);

% State vector X
% TVMAT must have increasing time grid points 
if isequal(Opt.Type,'ForwardInTime')
    X = tvmat(x,t);
else
    X = tvmat(flip(x),flip(t));
end

%% Construct Outputs
[U1,C,D] = evalt(U,C,D,X.Time);
Y = C*X+D*U1;

% Return X and Y in terms of input U time grid
% Discussed with Pete during 12/6/2019
[X,Y] = evalt(X,Y,U.Time);

%% LOCALderiv
function xdot = LOCALderiv(t,x,A,B,U)
[At,Bt,ut] = tvsubs(A,B,U,t);
xdot = At*x + Bt*ut;