function [t,X,Y] = simsfunc(mySystem,tSpan,varargin)
% SIMSFUNC Runs s-function system and calculate outputs.

%% Get Sizes
sizeInfo = feval(mySystem,[],[],[],0);
Nx = sizeInfo(1);  % Number of continuous states
Ny = sizeInfo(3);  % Number of outputs
Nu = sizeInfo(4);  % Number of inputs

%% Input Processing
narginchk(2,5);
nin = nargin;

% Defaults
x0        = zeros(Nx,1);
Opt       = tvodeOptions;
Input     = tvmat(zeros(Nu,1,length(tSpan)),tSpan);

switch nin
    case 3
        Input = varargin{1};
    case 4
        Input = varargin{1};
        x0    = varargin{2};
    case 5
        Input = varargin{1};
        x0    = varargin{2};
        Opt   = varargin{3};
end

% Read Opts
OdeSolver = Opt.OdeSolver;
OdeOpt    = Opt.OdeOptions;

%% Define S-funciton flags
derivativesFlag = 1;
ouputsFlag = 3;

%% Simulate ODE
odefh = str2func(OdeSolver);
[t,x] = odefh(@NESTED_SDRIVER,tSpan,x0,OdeOpt);
X = tvmat(x,t);

%% NESTED function used by ODE solver to 'drive' the s-function
    function xdot = NESTED_SDRIVER(t,x)
        u = tvsubs(Input,t);
        xdot = feval(mySystem,t,x,u(:),derivativesFlag);
    end

%% Process Outputs
Nt = length(t);
y = zeros(Ny,Nt);
for i = 1:Nt
    y(:,i) = feval(mySystem,t(i),x(i,:)',tvsubs(Input,t(i)),ouputsFlag);
end
Y = tvmat(y,t);
end