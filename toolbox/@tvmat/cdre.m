function [P,G,Pdot,sol] = cdre(A,B,Q,varargin)
% cdre  Solve Continuous-time Differential Riccati Equation.
%
% [P,G,Pdot,sol] = cdre(A,B,Q,R,S,E,F,Tspan) integrates the
% continuous-time Riccati differential equation:
%    Pdot + A'PE + E'PA - (E'PB + S) inv(R) (B'PE + S') + Q = 0
%    P(Tspan(1)) = F
% When omitted, R, S, and E are set to the default values R=I, S=0, and
% E=I. Besides the solution P, cdre also returns the gain matrix:
%    G = inv(R) (E'PB + S)'
% The CDRE is integrated on the interval [Tspan(1) Tspan(end)]. If Tspan
% is a 2-by-1 vector then outputs (P,G,Pdot) are TVMATs specified on the
% solution times returned by ODE45. If Tspan is an Nt-by-1 vector (with
% Nt>2)  then the output TVMATs are specified on Tspan. Finally, sol is
% the ODE45 solutions structure. This can be used to evaluate the CDRE
% solution at any time t using DEVAL:
%    [P,Pdot] = deval(sol,t);
%
% [P,G,Pdot,sol] = cdre(A,B,Q,...,Tspan,Opt) allows options for the ODE
% solver to be specified.  Opt is a structure with 2 fields. OdeSolver and
% OdeOptions. See ODESET for more help.
%
% (Optional) Options if specified then must be of type tvodeOptions which
% can be provided and customized using
% cdreOpt = tvodeOptions;
% cdreOpt.OdeOptions = odeset('RelTol',1e-5);
%
% If Options are not specified then cdre considers following default options.
% OdeSolver = 'ode45';
% OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);

% XXX This function should have hooks to allow (A,B,Q,...) to be function
% handles (since that is all that is required to do the integration)

%% Input Processing
%  [P,G,Pdot,sol] = cdre(A,B,Q,R,S,E,F,Tspan,opt)
nin = nargin;
narginchk(5,9);
[R,S,E,F,Tspan,opt] = LOCALinputprocessing(nin,varargin);
A = tvmat(A);
B = tvmat(B);
Q = tvmat(Q);
InterpMethod = A.InterpolationMethod;

%% Set Default Values
Nx = size(A,1);
if isempty(F)
    F = zeros(Nx);
    F = F(:);
end

% PARSE opt to OdeSolver and OdeOptions
if isempty(opt) % If isempty(opt) then specify default values for OdeSolver and
    opt = tvodeOptions('OdeOptions',odeset('RelTol',1e-5,'AbsTol',1e-8));
end
% XXX - The code below overwrites any Event specified by the user.
OdeOpt = odeset(opt.OdeOptions,'Events',@LOCALevents);
OdeSolverStr = opt.OdeSolver;
OdeSolver = str2func(OdeSolverStr);

%% Create function handles for ODE45
[Afh,Bfh,Qfh,Rfh,Sfh,Efh] = tv2fh(A,B,Q,R,S,E);
odefh = @(t,P) LOCALPdot(t,P,Afh,Bfh,Qfh,Rfh,Sfh,Efh);

%% Solve Riccati Differential Equation
% Note: Setting the time span as [T, 0] indicates that OdeSolver should
% integrate backwards from the boundary condition P(T)=F.
% warning off;
%[t,P] = ode45( odefh,tspan,F,odeopts);
warning('off',['MATLAB:' OdeSolverStr ':IntegrationTolNotMet']);
sol = OdeSolver(odefh,Tspan,F,OdeOpt);
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
    Nt = numel(t);
end

% Compute P and Pdot
nout = nargout;
if nout<=2
    P = deval(sol,t);
    P = tvmat( reshape(P,[Nx Nx Nt]) , t, InterpMethod);
else
    % Two options for computing Pdot:
    %  A) DEVAL computes Pdot as "the first derivative of the polynomial
    %        interpolating the solution"
    %  B) Directly call the ODEFH to compute Pdot given P
    % Option B seems like it should be more accurate (but perhaps
    % more computationally costly?)
    
    [P,Pdot] = deval(sol,t);
    P = tvmat( reshape(P,[Nx Nx Nt]) , t , InterpMethod);
    Pdot = tvmat( reshape(Pdot,[Nx Nx Nt]) , t , InterpMethod);
    
    %     Pdot = zeros(Nx,Nx,Nt);
    %     for i=1:Nt
    %         Pdot(:,:,i) = reshape( odefh(t(i),P(:,:,i)), [Nx Nx]);
    %     end
    %     Pdot = tvmat(Pdot,t);
    
    % Enforce Symmetry
    Pdot = (Pdot+Pdot')/2;
end

% Enforce Symmetry
P = (P+P')/2;

% Compute Gain
if nout>=2
    B = evalt(B,t);
    
    if ~isempty(E)
        E = evalt(E,t);
        EPB = E'*P*B;
    else
        EPB = P*B;
    end
    
    if ~isempty(S)
        S = evalt(S,t);
        G = (EPB + S)';
    else
        G = EPB';
    end
    
    if ~isempty(R)
        R = evalt(R,t);
        G = R\G;
    end
end

%% LOCAL Function
% This can be improved by:
%  1) Accounting for the symmetry in P;
%  2) Possibly doing some initial transformation on (A,B,Q,S,R)
%     to put the data in a form for more efficient integration.
%     [Perhaps look at work by Laub, others for AREs]
%  3) Connection to integrating Hamiltonian?  Look at the PLTV literature.
function Pdot = LOCALPdot(t,P,Afh,Bfh,Qfh,Rfh,Sfh,Efh)

% Get RDE matrices
A = Afh(t);
B = Bfh(t);
Q = Qfh(t);

% Convert P from column to matrix
Nx = size(A,1);
P = reshape(P,[Nx, Nx]);

% Compute Pdot (Defaults are R=I, S=0, and E=I)
if isempty(Efh)
    EPA = P*A;
    EPBS = P*B;
else
    E = Efh(t);
    EPA = E'*P*A;
    EPBS = E'*P*B;
end
if ~isempty(Sfh)
    S = Sfh(t);
    EPBS = EPBS + S;
end
if isempty(Rfh)
    Pdot = -EPA -EPA' - Q + EPBS*EPBS';
else
    R = Rfh(t);
    Pdot = -EPA -EPA' - Q + EPBS/R*EPBS';
end

% Enforce Symmetry
Pdot = (Pdot+Pdot')/2;

% Convert Pdot from matrix back to column
Pdot = Pdot(:);

%% LOCAL Function
function [value,isterminal,direction] = LOCALevents(~,P)

if any(isnan(P(:)))
    value = 0;
else
    value = 1;
end
isterminal = 1;
direction = 0;

%% LOCAL Function
function [R,S,E,F,Tspan,Opt] = LOCALinputprocessing(nin,vin)

R=[]; S=[]; E=[]; Opt=[];
if nin==5
    % cdre(A,B,Q,F,Tspan)
    [F,Tspan] = deal(vin{:});
elseif nin==6
    if isa(vin{end},'struct')
        % cdre(A,B,Q,F,Tspan,opt)
        [F,Tspan,Opt] = deal(vin{:});
    else
        % cdre(A,B,Q,R,F,Tspan)
        [R,F,Tspan] = deal(vin{:});
    end
elseif nin==7
    if isa(vin{end},'struct')
        % cdre(A,B,Q,R,F,Tspan,opt)
        [R,F,Tspan,Opt] = deal(vin{:});
    else
        % cdre(A,B,Q,R,S,F,Tspan)
        [R,S,F,Tspan] = deal(vin{:});
    end
elseif nin==8
    if isa(vin{end},'struct')
        % cdre(A,B,Q,R,S,F,Tspan,opt)
        [R,S,F,Tspan,Opt] = deal(vin{:});
    else
        % cdre(A,B,Q,R,S,E,F,Tspan)
        [R,S,E,F,Tspan] = deal(vin{:});
    end
elseif nin==9
    % cdre(A,B,Q,R,S,E,F,Tspan,opt)
    [R,S,E,F,Tspan,Opt] = deal(vin{:});
end

% Lift to TVMATs
R = tvmat(R);
S = tvmat(S);
E = tvmat(E);
