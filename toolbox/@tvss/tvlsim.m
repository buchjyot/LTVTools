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
%        [Optional] x0 is initial condition, default to 0
%        [Optional] Opt is tvOdeOptions

%% Input Processing
narginchk(2,4);
nin = nargin;

Nx = order(G);
Opt = [];x0 = [];
switch nin
    case 3
        if isa(varargin{1},'tvodeOptions')
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
    % Assume zero initial conditions if user did not provided any i.e.
    % tvlsim(G,U,[],Opt)
    x0 = zeros(Nx,1);
end
if isempty(Opt)
    Opt = tvodeOptions('OdeOptions',odeset('RelTol',1e-5,'AbsTol',1e-8));
end
OdeSolverFh = str2func(Opt.OdeSolver);
OdeOpt = Opt.OdeOptions;

% Make sure the system (G) is defined on the given input (U) horizon
[A,B,C,D] = ssdata(G);
ltvutil.verifyFH(G,U);
[GT0,GTf] = getHorizon(G);
[UT0,UTf] = getHorizon(U);
if any(UT0>GTf || UTf>GTf || UT0<GT0 || UTf<GT0)
    error('Input U must be defined on the subset of the system horizon.');
end

%% Simulate System
% Solve ODE on [GT0,GTf]
odefh = @(t,x) LOCALderiv(t,x,A,B,U);
[t,x] = OdeSolverFh(odefh,[GT0 GTf],x0,OdeOpt);

% State vector X
X = tvmat(x,t);

%% Construct Outputs
[U,C,D] = evalt(U,C,D,t);
Y = C*X+D*U;

%% LOCALderiv
function xdot = LOCALderiv(t,x,A,B,U)
[At,Bt,ut] = tvsubs(A,B,U,t);
xdot = At*x + Bt*ut;