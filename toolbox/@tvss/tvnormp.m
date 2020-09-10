function [glb,dwc,info] = tvnormp(G,NE,U1,x0,Opt)
%% TVNORMP
% Runs power iterations to estimate the best lower bound
%
% Input argument are as follows:
% G : LTV System
% NE: Number of Euclidean outputs
% d1: starting input for power iterations
% Opt: poweritOptions
%
% Outputs are power iteration lower bound, worst-case disturbance and info
% structure.

% Input processing
narginchk(1,5);
nin = nargin;
switch nin
    case 1
        NE = 0;U1 = [];x0 = [];Opt = poweritOptions;
    case 2
        U1 = [];x0 = [];Opt = poweritOptions;
    case 3
        x0 = [];Opt = poweritOptions;
    case 4
        Opt = poweritOptions;
end

% Precompute Adjoint System
[A,B,C,D] = ssdata(G);
Nu = size(B,2);
[NY,Nx] = size(C);
NL2 = NY-NE;
Ga = tvss(-A',-C(1:NL2,:)',B',D(1:NL2,:)');
[T0,Tf] = getHorizon(G);
CTf = tvsubs(C(NL2+1:end,:),Tf);

% Initial condition and input
if isempty(x0)
    x0 = zeros(Nx,1);
end
if isempty(U1)
    U1 = tvmat(randn(Nu,1,20),linspace(T0,Tf,20));
end

% Options
MaxIter    = Opt.MaxIter;
StopTol    = Opt.StopTol;
OdeSolver  = Opt.OdeSolver;
OdeOptions = Opt.OdeOptions;
Display    = isequal(Opt.Display,'on');
odeOpt     = tvodeOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions);

% Initial Input
U1 = U1/tvnorm(U1);

% Display
if Display
    fprintf(' ### Starting Power Iterations:\n');
end

% Call the solver
[glb,dwc,info] = ltvutil.tvpowerit(G,Ga,U1,CTf,T0,Tf,NL2,x0,MaxIter,StopTol,Display,odeOpt);
end