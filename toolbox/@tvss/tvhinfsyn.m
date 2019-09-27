function [K,CL,g,info] = tvhinfsyn(G,Ny,Nu,varargin)
%% TVHINFSYN  H-infinity measurement-feedback synthesis on finite horizon.
%
%   [K,CL,GAM] = TVHINFSYN(G,NMEAS,NCON,NEUCLIDEAN) calculates the
%   Hoo-optimal control law u = K y for the LTV plant G with state-space
%   equations:
%
%       dx =  A(t) x +  B1(t) w +  B2(t) u
%      eL2 = C1(t) x + D11(t) w + D12(t) u
%       eE = CE(t) x + DE(t) u
%        y = C2(t) x + D21(t) w + D22(t) u
%
%   NMEAS and NCON specify the numbers of measurements y and controls u (u
%   and y must be the last inputs and outputs of G). NEUCLIDEAN specify
%   Euclidean number of penalty on outputs. The controller K minimizes the
%   finite horizon induced-L2 and induced-L2-to-Euclidean norm of CL =
%   LFT(G,K), the closed-loop transfer function from disturbance signals w
%   to error signals z.
%
%   [K,CL,GAM] = TVHINFSYN(G,NMEAS,NCON,GAMTRY) calculates the H-infinity
%   controller for the performance level GAMTRY. When GAMTRY is feasible,
%   GAM is the actual closed-loop performance obtained with this controller
%   (GAM<=GAMTRY). Otherwise GAM=Inf and K and CL are set to [].
%
%   [K,CL,GAM] = TVHINFSYN(G,NMEAS,NCON,OPT) searches the range OPT.BOUNDS =
%   [GMIN,GMAX] for the best achievable performance GAM. TVHINFSYN returns a
%   controller K with performance
%      * GAM <= GMIN when GMIN is achievable
%      * GMIN < GAM <= GMAX when GMAX but not GMIN is achievable
%   If GMAX itself is not achievable, TVHINFSYN returns GAM=Inf and K=CL=[].
%   Use tvhinfsynOptions to create the option set OPT.
%
%   [K,...] = TVHINFSYN(G,NMEAS,NCON,OPT) specifies additional options.
%   Use tvhinfsynOptions to create the option set OPT.
%
%   [K,CL,GAM,INFO] = TVHINFSYN(G,NMEAS,NCON,...) also returns a structure
%   INFO with additional synthesis data. For Riccati-based synthesis,
%   this includes
%       gamma   Performance level used to compute K
%           X   Riccati solution Xoo for this performance level
%           Y   Riccati solution Yoo for this performance level
%       Ku,Kw   State-feedback gains
%       Lx,Lu   Observer gains
%        Preg   Regularized plant used to compute K
%          AS   All controllers with performance INFO.gamma are given by
%               K = lft(INFO.AS,Q) where Q is any stable transfer function
%               of size [NCON NMEAS] with peak gain less than INFO.gamma.
%   For LMI-based synthesis, INFO contains the best performance GAMMA and
%   the corresponding LMI solutions R and S.
%
%   The observer form of the controller K is
%       dxe = A xe + B1 we + B2 u + Lx e
%         u = Ku xe + Lu e
%        we = Kw xe
%   where we is an estimate of the worst-case perturbation and
%         e = y - C2 xe - D21 we - D22 u
%   plays the role of "innovation."
%
%   See also HINFSYNOPTIONS, HINFFI, HINFFC, HINFSTRUCT, H2SYN, LTRSYN.

%% List Assumptions

% Finite horizon generalized Hinf synthesis problem (Output Feedback):
%
% (A1) In oder that a stabilizing controller exists, we assume that the pair
% (A,B2) is stabilizable.
%
% (A2) Pair (A,C2) is detectable.
%
% (A3) D12'*D12 = eye(Nu) and D21*D21' = eye(Ny) have full rank for all
% times of interest.
%
% (A4) K is causal, linear, time-varying controller satisfying the objective
% ||CL|| < g on finite horizon [T0 Tf].

%% Existance of a Solution

% Hinf generalized regulator problem has a solution on a finite horizon if
% and only if,
%
% (1) The CDRE associated with full-information control problem has a
% solution on a finite horizon
%
% (2) The CDRE associated with Hinf estimation has a solution on the same
% finite horizon
%
% (3) Coupling condition is satisfied

%% Input Processing & Initialization

% Number of input and output arguments
nin = nargin;
nout = nargout;

% Minimum 3 and maximum 6 input arguments are allowed
narginchk(3,6);

% Initialize Input and Output Arguments
gTry = [];

% [gLow,gUpp] are bounds to start with for bisection
Opt = [];
NE = 0;
switch nin
    case 4
        if isa(varargin{1},'tvhinfsynOptions')
            Opt = varargin{1};
        elseif isa(varargin{1},'double')
            NE = varargin{1};
        end
    case 5
        NE = varargin{1};
        if isa(varargin{2},'tvhinfsynOptions')
            Opt = varargin{2};
        elseif isa(varargin{2},'double')
            gTry =  varargin{2};
        end
    case 6
        % Default Order
        NE = varargin{1};
        gTry = varargin{2};
        Opt = varargin{3};
end

% Use Default Options
if isempty(Opt)
    Opt = tvhinfsynOptions;
end

% Time Horizon
ltvutil.verifyFH(G);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% Check if NE is valid
NY = size(G,1);
if ~(NE >= 0) || ~(NE <= NY-Nu)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY-Nu.');
end

%% Call Engine
if nout~=4
    [K,CL,g] = ltvutil.tvhinfric(G,Ny,Nu,NE,gTry,Opt);
else
    [K,CL,g,info] = ltvutil.tvhinfric(G,Ny,Nu,NE,gTry,Opt);
end
end