function [W,Wdot,info] = cdle(A,B,W0,Tspan,varargin)
% cdle - Solve Continuous-time Diffrential matrix Lyapunov Equation
%
% [W,opt] = cdle(A,B) solves the following continuous-time Diffrential
% Matrix Lyapunov Equation
%
%           A*W + W*A' + B*B' = d(W)/dt
%           W(Tspan(1)) = F
%
% NOTE: cdle uses the default horizon as specified in the A and B matrices
% if Tspan is not specified by user

%% Input Processing
narginchk(2,5);
Nx = size(A,1);
ltvutil.verifyFH(A);
[T0,Tf] = getHorizon(A);

opt = [];
switch nargin
    case 2
        % Specify Tspan (Forward in Time by default)
        Tspan = [T0,Tf];
        % Specify Initial Condition as 0
        W0 = zeros(Nx);
    case 3
        % Specify Tspan (Forward in Time by default)
        Tspan = [T0,Tf];
    case 4
        opt = varargin{1};
end

% LIFT to tvmat
A = tvmat(A);
B = tvmat(B);
InterpMethod = A.InterpolationMethod;

%% Set Default Values
% PARSE opt to OdeSolver and OdeOptions
if isempty(opt) % If isempty(opt) then specify default values for OdeSolver and
    opt = tvodeOptions('OdeOptions',odeset('RelTol',1e-5,'AbsTol',1e-8));
end
% XXX - The code below overwrites any Event specified by the user.
OdeOpt = odeset(opt.OdeOptions,'Events',@LOCALevents);
OdeSolverStr = opt.OdeSolver;
OdeSolver = str2func(OdeSolverStr);

%% Create function handles for ODE solver
[Afh,Bfh] = tv2fh(A,B);
odefh = @(t,W) LOCALWdot(t,W,Afh,Bfh);

%% Solve Lyapunov Diffrential Equation
% Note: Setting the time span as [T, 0] indicates that OdeSolver should
% integrate backwards from the boundary condition P(T) = W0.
W0 = W0(:);
warning('off',['MATLAB:' OdeSolverStr ':IntegrationTolNotMet']);
sol = OdeSolver(odefh,Tspan,W0,OdeOpt);
warning('on',['MATLAB:' OdeSolverStr ':IntegrationTolNotMet']);

%% Construct Outputs
% Note: sol = ode45() returns a solution structure that can be
% evaluated with deval. The alternative syntax [t,P]=ode45() refines
% the solution time grid by a factor of 4.  In other words,
%   t(1)=sol.x(1), t(5)=sol.x(2), t(9)=sol.x(3).
% The documentation states this default refinement is used because
% ODE45 can take large time steps.

% Time Vector
if numel(Tspan)==2
    if isempty(OdeOpt.Refine)
        Nr = 4;    % This is the default ODESET option
    else
        Nr = OdeOpt.Refine;
    end
    fac = ((1:Nr)-1)/Nr;
    tSol = sort(sol.x);
    Nsol = numel(tSol);
    tDiff = diff(tSol);
    t = zeros(Nr,Nsol-1);
    t(1,:) = tSol(1:end-1);
    for i=2:Nr
        t(i,:) = tSol(1:end-1) + fac(i)*tDiff;
    end
    t = [t(:); tSol(end)];
    Nt = numel(t);
else
    t = sort(Tspan);
end

% Compute W and Wdot
nout = nargout;
if isequal(nout,1)
    W = deval(sol,t);
    W = tvmat( reshape(W,[Nx Nx Nt]) , t , InterpMethod);
else
    % Two options for computing Pdot:
    %  A) DEVAL computes Pdot as "the first derivative of the polynomial
    %        interpolating the solution"
    %  B) Directly call the ODEFH to compute Pdot given P
    % Option B seems like it should be more accurate (but perhaps
    % more computationally costly?)
    
    [W,Wdot] = deval(sol,t);
    W = tvmat( reshape(W,[Nx Nx Nt]) , t , InterpMethod);
    Wdot = tvmat( reshape(Wdot,[Nx Nx Nt]) , t , InterpMethod);
    Wdot = (Wdot+Wdot')/2;
end
info.sol = sol;

% Brute Force Symmetry to avoid numerical issues
W = (W+W')/2;
end

%% LOCAL Function
% This can be improved by:
%  1) Accounting for the symmetry in W;
function Wdot = LOCALWdot(t,W,Afh,Bfh)

% Get Lyapunov state matrices
A = Afh(t);
B = Bfh(t);

% Convert W from column to matrix
Nx = size(A,1);
W = reshape(W,[Nx, Nx]);

% Compute Wdot (Lyapunov Diffrential Equation)
Wdot = A*W + W*A' + B*B';
Wdot = (Wdot+Wdot')/2;

% Convert Pdot from matrix back to column
Wdot = Wdot(:);
end

%% LOCAL Function
function [value,isterminal,direction] = LOCALevents(~,W)

if any(isnan(W(:)))
    value = 0;
else
    value = 1;
end
isterminal = 1;
direction = 0;
end