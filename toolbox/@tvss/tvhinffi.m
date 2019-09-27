function [K,CL,g,info] = tvhinffi(G,Nu,varargin)
%% TVHINFFI  Full-information H-infinity Synthesis on a Finite Horizon.
%
%   Full-information synthesis assumes the controller has access to both
%   the state vector x and disturbance w. Use HINFSYN for the more general
%   output-feedback case when only output measurements are available.
%
%   [K,CL,GAM] = TVHINFFI(G,NCON) calculates the Hoo-optimal control law
%   u = K [x;w] for the LTV plant G with state-space equations:
%
%       dx =   A(t) x +  B1(t) w +  B2(t) u
%        z =  C1(t) x + D11(t) w + D12(t) u
%
%   NCON specifies the numbers of controls u (must be the last inputs of G).
%   The gain matrix K minimizes the H-infinity norm of the closed-loop
%   transfer function from disturbance signals w to error signals z.
%
%   [K,CL,GAM] = TVHINFFI(G,NCON,GAMTRY) calculates a gain matrix K that
%   achieves the closed-loop performance level GAMTRY. If GAMTRY is not
%   achievable, TVHINFFI returns GAM=Inf and K=CL=[].
%
%   [K,CL,GAM] = TVHINFFI(G,NCON,OPT) searches the range [GMIN,GMAX]
%   for the best achievable performance GAM. TVHINFFI returns a gain matrix
%   K with performance
%      * GAM <= GMIN when GMIN is achievable
%      * GMIN < GAM <= GMAX when GMAX but not GMIN is achievable
%   If GMAX itself is not achievable, TVHINFFI returns GAM=Inf and K=CL=[].
%  [GMIN,GMAX] are specified using tvhinfsynOptions Bounds property.
%
%   [K,...] = TVHINFFI(G,NCON,...,OPT) specifies additional options.
%   Use tvhinfsynOptions to create the option set OPT.
%
%   [K,CL,GAM,INFO] = TVHINFFI(G,NCON,...) also returns a structure INFO
%   with the following synthesis data:
%        gamma   Performance level used to compute K
%            P   Riccati solution Poo for this performance level
%         Greg   Regularized plant used to compute K
%
%   Reference: Green, M., & Limebeer, D. J. (2012). Linear robust control.
%   Courier Corporation. Chapter 6
%
%   See also HINFSYNOPTIONS, HINFFC, HINFSYN.

%% List Assumptions

% Finite horizon full-information Hinf synthesis problem:
%
% (1)In oder that a stabilizing controller exists, we assume that the pair
% (A,B2) is stabilizable.
%
% (2) Pair (A,C1) has no unobservable mode on the imaginary axis.
%
% (3) K is causal, linear, time-varying controller satisfying the objective
% ||CL|| < g on finite horizon [T0 Tf].

%% Input Processing & Initialization

% Number of input and output arguments
nin = nargin;
nout = nargout;

% Minimum 2 and maximum 5 input arguments are allowed
narginchk(2,5);

% Initialize Input and Output Arguments
gTry = [];

% [gLow,gUpp] are bounds to start with for bisection
Opt = [];
NE = 0;
switch nin
    case 3
        if isa(varargin{1},'tvhinfsynOptions')
            Opt = varargin{1};
        elseif isa(varargin{1},'double')
            NE = varargin{1};
        end
    case 4
        NE = varargin{1};
        if isa(varargin{2},'tvhinfsynOptions')
            Opt = varargin{2};
        elseif isa(varargin{2},'double')
            gTry =  varargin{2};
        end
    case 5
        % Default Order
        NE = varargin{1};
        gTry = varargin{2};
        Opt = varargin{3};
end

% Use Default Options
if isempty(Opt)
    Opt = tvhinfsynOptions;
end

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% Check if NE is valid
NY = size(G,1);
if ~(NE >= 0) || ~(NE <= NY-Nu)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY-Nu.');
end

% Verify Finite Horizon System
ltvutil.verifyFH(G);

%% Call Engine
if nout~=4
    [K,CL,g] = ltvutil.tvhinfKFI(G,Nu,NE,gTry,Opt);
else
    [K,CL,g,info] = ltvutil.tvhinfKFI(G,Nu,NE,gTry,Opt);
end
end