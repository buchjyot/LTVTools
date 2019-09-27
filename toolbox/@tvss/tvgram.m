function [W,Wdot,info] = tvgram(sys,b,varargin)
% TVGRAM Controllability and observability gramians for LTV Systems on
% Finite Horizon
%
%   Wc = TVGRAM(SYS,'c') computes the controllability gramian of the
%   state-space model SYS (see TVSS) over the horizon as default time
%   scale.
%
%   Wo = TVGRAM(SYS,'o') computes its observability gramian.
%
%   In both cases, the state-space model SYS should be stable. The gramians
%   are computed by solving continuous-time matrix Lyapunov Diffrential
%   Equation:
%
%     A*Wc + Wc*A' + B*B' = d(Wc)/dt   [Wc is controllability gramian]
%     A'*Wo + Wo*A + C'*C = -d(Wo)/dt   [Wo is observability gramian]
%
%     for continuous-time Linear Time-Varying systems
%         dx/dt = A(t) x(t) + B(t) u(t),   y = C(t) x(t) + D(t) u(t)
%
% For a discrete-time, time-varying model SYS:  XXX Not Implemented
%
% Example:
%
%   Time = linspace(0,5,10)';
%   AData = -5+0.1*Time.^2;
%   A = tvmat(AData,Time);
%   B = 1; C = 1; D=0;
%   G = tvss(A,B,C,D);
%   Wc = tvgram(G,'c');
%   % Options can be specified using tvodeOptions
%   tvgramopt = tvodeOptions;
%   tvgramopt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
%   Wo = tvgram(G,'o',tvgramopt);

% Input Processing
narginchk(2,3);
if ~isa(sys,'tvss') && ~isa(sys,'tvuss')
    error('Input System must be Time Varying State Space');
end
if ~ischar(b)
    error('Second Input must be char');
end
if nargin>2
    cdleopt = varargin{1};
else
    cdleopt = [];
end

% Extract ssdata from G
[A,B,C,D] = ssdata(sys);
if any(D.Data(:))
    error('Feedthrough matrix D must be zero.');
end

% Time Horizon
ltvutil.verifyFH(sys);
[T0,Tf] = getHorizon(sys);

% Initialize W(Tf) as F being all zeros
Nx = size(A,1);
F = zeros(Nx);

% Compute the gramians based on Theorem 3.3.1 from Green and Limebeer,
% Linear Robust Control
switch lower(b)
    case 'c'
        % Solve cdle A*Wc + Wc*A' + BB' = d(Wc)/dt
        [W,Wdot,info] = cdle(A,B,F,[T0 Tf],cdleopt);
    case 'o'
        % Solve cdle A'*Wo + Wo*A + C'C = -d(Wo)/dt
        [W,Wdot,info] = cdle(-A',1i*C',F,[Tf T0],cdleopt);
        Wdot = -1*Wdot;W = -W;info.sol.y = -1*info.sol.y;
    otherwise
        error('Second Input must be either ''c'' for controllability gramian or ''o'' for observability gramian.')
end
end