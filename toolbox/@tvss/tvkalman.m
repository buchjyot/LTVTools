function [Kest,L,P,Pdot,sol] = tvkalman(SYS,Qn,Rn,varargin)
% TVKALMAN Time-Varying Kalman State Estimator
%
%   [KEST,L,P] = TVKALMAN(SYS,QN,RN,NN) designs a Kalman estimator KEST for
%   continuous-time plants
%      .
%      x = A(t)x(t) + B(t)u(t) + G(t)*w        {State equation}
%      y = C(t)x(t) + D(t)u(t) + H(t)*w + v    {Measurements}
%
%   with known inputs u(t), process noise w, and measurement noise v,
%   KEST uses [u(t);y(t)] to generate optimal estimates y_e(t),x_e(t) of
%   y(t),x(t) by:
%       .
%      x_e  = A(t)x_e(t) + B(t)u(t) + L(t) (y - C(t)x_e(t) - D(t)u(t))
%
%      |y_e(t)| = |C(t)| x_e(t) + |D(t)| u(t)
%      |x_e(t)| = |I|    x_e(t) + |0|    u(t)
%
%   TVKALMAN takes the state-space model SYS = TVSS(A,[B G],C,[D H]) and
%   the constant covariance matrices:
%
%      QN = E{ww'},     RN = E{vv'},     NN = E{wv'}.
%
%   The row size of QN specifies the length of w and NN is set to 0 when
%   omitted. KALMAN returns the estimator gain L and the steady-state error
%   covariance P (solution of the associated Riccati equation).
%
%   For discrete-time plants, XXX TVKALMAN is not implemented yet

%% Input Processing
narginchk(5,7);

% Convert to strings to chars as strcmp does not work on cell arrays that
% have strings in them.
varargin = controllib.internal.util.hString2Char(varargin);

% tvkalman process only one system at a time
if ndims(SYS)>2 %#ok<ISMAT>
    error(message('Control:general:RequiresSingleModel','tvkalman'));
end

% Get the inputs
nin = nargin;
vin = varargin;

% Extract Plant Data
[A,BB,C,DD] = ssdata(SYS);

% A is a square matrix and DD is an augmented matrix [D H]
Nx = size(A,1); %#ok<*NASGU>
[No,NupNw] = size(DD);

% Process varargin
E = [];  % Descriptor Systems not currently handled.
Nn = [];
Opt = [];
if nin==5
    % tvkalman(sys,Q,R,F,Tspan)
    [F,Tspan] = deal(vin{:});
elseif nin==6
    if isa(vin{end},'tvodeOptions')
        % tvkalman(sys,Q,R,F,Tspan,opt)
        [F,Tspan,Opt] = deal(vin{:});
    else
        % tvkalman(sys,Q,R,N,F,Tspan)
        [Nn,F,Tspan] = deal(vin{:});
    end
else
    [Nn,F,Tspan,Opt] = deal(varargin{:});
end
% NOTE: If F is empty then it will be initialized to zeros by cdre function
% F in tvkalman represents the initial condition error covariance

% Validate Qn,Rn,Nn
[Nw,Nw2] = size(Qn);
[Ny,Ny2] = size(Rn);
Nu = NupNw-Nw;
if Nw~=Nw2 || Nw>NupNw || ~isreal(Qn)
    error(message('Control:design:kalman2','QN',NupNw))
elseif Ny~=Ny2 || Ny>No || ~isreal(Rn)
    error(message('Control:design:kalman2','RN',No))
elseif Ny==0
    error(message('Control:design:kalman1'))
elseif nin<5 && Ny~=No
    error(message('Control:design:kalman3'))
elseif (~isequal(size(Nn),[Nw Ny]) || ~isreal(Nn)) && ~isempty(Nn)
    error(message('Control:design:kalman4'))
elseif norm(Qn'-Qn,1) > 100*eps*norm(Qn,1)
    warning(message('Control:design:MakeSymmetric','tvkalman(SYS,Qn,Rn,Nn)','Qn','Qn','Qn'))
elseif norm(Rn'-Rn,1) > 100*eps*norm(Rn,1)
    warning(message('Control:design:MakeSymmetric','tvkalman(SYS,Qn,Rn,Nn)','Rn','Rn','Rn'))
end

% Extract B,G,D,H
B = BB(:,1:Nu);
D = DD(:,1:Nu);
G = BB(:,Nu+1:end);
H = DD(:,Nu+1:end);

% XXX - Update this (Force H to be zero)
H = 0;

%% Update Noise Covariances
Qnbar = G*Qn*G';
if ~isempty(Nn)
    Rnbar = Rn + H*Nn + Nn'*H' + H*Qn*H';
    Nnbar = G*(Qn*H' + Nn);
elseif ~isequal(H,0)
    Rnbar = Rn + H*Qn*H';
    Nnbar = G*(Qn*H');
else
    Rnbar = Rn;
    Nnbar = Nn;
end

%% Enforce Duality
% To enforce duality and to use CDRE transpose A and C
At = A';
Ct = C';

% NOTE: CDRE function solves a problem like
%       Pdot = -A'*P - P*A - Q + EPBS*EPBS';
%
%       Riccati Diffrential Equation for LQR
%       -Pdot = A'*P + P*A - P*B*inv(R)*B'*P + G*Q*G'
%       Pdot = -A'*P - P*A + P*B*inv(R)*B'*P - G*Q*G'  {Same as CDRE}
%
%       Riccati Diffrential Equation for Kalman Filter
%       Pdot = A*P + P*A' - P*C'*inv(R)*C*P + G*Q*G'
%       Pdot = -(-A)*P - P*(-A') + P*C'*inv(-R)*C*P - (-G*Q*G') {Same as CDRE}
At = -At;
Rnbar = -Rnbar;
Qnbar = -Qnbar;
Nnbar = -Nnbar;

%% Keep the Time Vector as it is
% CDRE Convention is P(T0)=F so as the TVKALMAN convention
% NOTE: F is P(t0) i.e. initial condition error covariance
Tspan = sort(Tspan(:));

%% Solve the CDRE associated with the time-varying Kalman Filter problem.
nout=nargout;
if nout<=2
    [P,L] = cdre(At,Ct,Qnbar,Rnbar,Nnbar,E,F,Tspan,Opt);
else
    [P,L,Pdot,sol] = cdre(At,Ct,Qnbar,Rnbar,Nnbar,E,F,Tspan,Opt);
end

% Negate the gain as an effect of (-R)
% Transpose the output to honour duality
L = -L';

%% Create Estimator TVSS
%   KEST uses [u(t);y(t)] to generate optimal estimates y_e(t),x_e(t) of
%   y(t),x(t) by:
%       .
%      x_e  = A(t)x_e(t) + B(t)u(t) + L(t) (y - C(t)x_e(t) - D(t)u(t))
%
%      |y_e(t)| = |C(t)| x_e(t) + |D(t)| u(t)
%      |x_e(t)| = |I|    x_e(t) + |0|    u(t)
L = evalt(L,SYS.Time);
Ae = A-L*C;
Be = [B-L*D L];
Ce = [C;eye(Nx)];
De = [D zeros(Ny);zeros(Nx,Nu+Ny)];
Kest = tvss(Ae,Be,Ce,De);