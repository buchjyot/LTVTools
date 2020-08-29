function [A,B,C,D] = linsfunc(mySystem,varargin)
% LINSFUNC  Linearize s-function system about an operating point.
%   [A,B,C,D]=LINSFUNC(SYS)  Linearize system about origin with no input.
%
%   [A,B,C,D]=LINSFUNC(SYS,X)  Linearize system about X with no input.
%
%   [A,B,C,D]=LINSFUNC(SYS,X,U)  Linearize system about operating point
%   (X,U).
%
%   [A,B,C,D]=LINSFUNC(SYS,X,U,PARAM)  Linearize system about operating
%   point (X,U) using parameters specified in PARAM.  PARAM(1) = delta
%   (perturbation size).  PARAM(2) = t (time of linearization).
%
% See also LINMOD.

%% Define constants
DEFAULT_DELTA = 1e-5;

%% Get system info
[sizeInfo,x0] = feval(mySystem,[],[],[],0);
nStates  = sizeInfo(1);  % Number of continuous states
nOutputs = sizeInfo(3);  % Number of outputs
nInputs  = sizeInfo(4);  % Number of inputs

%% See if user specified additional arguments
if nargin > 1
    x = varargin{1}(:);
else
    x = x0;
end
if nargin > 2
    u = varargin{2}(:);
else
    u = zeros(nInputs,1);
end
if nargin > 3
    param = varargin{3};
    if length(param) >= 2 && isnumeric(param)
        delta = param(1);
        t = param(2);
    else
        error('Parameter argument must be a 3-by-1 double array')
    end
else
    t = 0;
    delta = DEFAULT_DELTA;
end

%% Check integrity of parameters
if ~isscalar(t) || t < 0
    error('Time must be a positive scalar.')
end
if ~isscalar(delta) || delta < 0
    error('Perturbation size must be a positive scalar')
end
if length(x) ~= nStates
    error('Operating point has wrong number of states.')
end
if length(u) ~= nInputs
    error('Operating point has wrong number of inputs.')
end

%% Get function values at this operating point
f = feval(mySystem,t,x,u,1);
g = feval(mySystem,t,x,u,3);

%% Preallocate memory
A = zeros(nStates,nStates);
C = zeros(nOutputs,nStates);

%% Perturb each state in x
for iState = 1:nStates
    
    xPert = x;
    xPert(iState) = x(iState) + delta;
    %xPert(iState) =  delta + 1e-3*delta*abs(x(iState));
    
    dfdxi = ( feval(mySystem,t,xPert,u,1) - f ) / delta;
    A(1:nStates,iState) = dfdxi;
    
    dgdxi = ( feval(mySystem,t,xPert,u,3) - g ) / delta;
    C(1:nOutputs,iState) = dgdxi;
end

%% Preallocate memory
B = zeros(nStates,nInputs);
D = zeros(nOutputs,nInputs);

%% Perturb each input in u
for jInput = 1:nInputs
    
    uPert = u;
    uPert(jInput) = u(jInput) + delta;
    %uPert(jInput) = delta + 1e-3*delta*abs(u(jInput)) ;
    
    dfdui = ( feval(mySystem,t,x,uPert,1) - f ) / delta;
    B(1:nStates,jInput) = dfdui;
    
    dgdui = ( feval(mySystem,t,x,uPert,3) - g ) / delta;
    D(1:nOutputs,jInput) = dgdui;
end
