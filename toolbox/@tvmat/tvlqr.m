function [K,P,Pdot,sol] = tvlqr(A,B,Q,R,varargin)
% tvlqr  Linear-quadratic regulator design for time-varying systems.
%
%   [K,P,Pdot] = tvlqr(SYS,Q,R,N,F,Tspan) calculates the optimal 
%   time-varying gain matrix K such that:
%
%     * For a continuous-time, time-varying model SYS, the state-feedback
%       law u(t) = -K(t)x(t)  minimizes the cost function
%             J = x(T)'F x(T) + Integral_T0^Tf {x'Qx + u'Ru + 2*x'Nu} dt
%             subject to the system dynamics  dx/dt = Ax + Bu.
%       where Tspan = [T0 Tf].
%
%     % For a discrete-time, time-varying model SYS:  XXX Not Implemented
%
%     The matrix N is set to zero when omitted. Also returned are (P,Pdot)
%     the solution and its derivative of the associated continuous-time 
%     differential Riccati equation (CDRE). If Tspan is a 2-by-1 vector 
%     then outputs (K,P,Pdot) are TVMATs specified on the solution times 
%     returned by ODE45.  If Tspan is an Nt-by-1 vector (with Nt>2) then 
%     the output TVMATs are specified on Tspan.
%
%     [K,P,Pdot,sol] = tvlqr(SYS,Q,...,Tspan,Opt) allows options for the
%     ODE45 solver to be specified. See ODESET for more help. This syntax
%     also returns the ODE45 solutions structure sol. This can be used to 
%     evaluate the CDRE solution at any time t using DEVAL:
%        [P,Pdot] = deval(sol,t);
%
%     [K,P,Pdot,sol] = tvlqr(A,B,Q,R,N,F,Tspan,Opt) is an equivalent syntax
%     for time-varying models with dynamics dx/dt = Ax + Bu.


%% Input Processing

% Convert to strings to chars as strcmp does not work on cell arrays that
% have strings in them.
varargin = controllib.internal.util.hString2Char(varargin); 

narginchk(6,8);
nin = nargin;
vin = varargin;

E = [];  % Descriptor Systems not currently handled.
N = [];
Opt = [];
if nin==6
    % tvlqr(A,B,Q,R,F,Tspan)
    [F,Tspan] = deal(vin{:});
elseif nin==7
    if isa(vin{end},'struct')
        % tvlqr(A,B,Q,R,F,Tspan,opt)
        [F,Tspan,Opt] = deal(vin{:});
    else
        % tvlqr(A,B,Q,R,N,F,Tspan)
        [N,F,Tspan] = deal(vin{:});
    end
else
    [N,F,Tspan,Opt] = deal(varargin{:});
end

%% Flip Time Vector
% CDRE Convention is P(T0)=F while TVLQR convention is P(Tf)=F
Tspan = sort(Tspan(:));
Tspan = flip(Tspan);

%% Solve the CDRE associated with the time-varying LQR problem.
nout=nargout;
if nout<=2
    [P,K] = cdre(A,B,Q,R,N,E,F,Tspan,Opt);
else
    [P,K,Pdot,sol] = cdre(A,B,Q,R,N,E,F,Tspan,Opt);
end
