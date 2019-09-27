function [Klqg] = tvlqg(SYS,Q,R,Qn,Rn,varargin)
% TVLQG  Synthesis of Time-Varying LQG regulators and servo-controllers.
% XXX - update the following 
% 
%   KLQG = LQG(SYS,QXU,QWV) computes an optimal linear-quadratic Gaussian
%   (LQG) regulator KLQG given a state-space model SYS of the plant and
%   weighting matrices QXU and QWV.  The dynamic regulator KLQG uses the
%   measurements y to generate a control signal u that regulates y around
%   the zero value. Use positive feedback to connect this regulator to the
%   plant output y.
%
%                                     w |               | v
%                                       |   .-------.   |
%                       .--------.      '-->|       |   V
%               .------>|  KLQG  |--------->|  SYS  |---O--.-----> y
%               |       '--------'  u       |       |      |
%               |                           '-------'      |
%               |                                          |
%               '------------------------------------------'
%
%   The LQG regulator minimizes the cost function
%
%         J(u) = Integral [x',u'] * QXU * [x;u] dt
%
%   subject to the plant equations
%
%         dx/dt = Ax + Bu + w
%             y = Cx + Du + v
%
%   where the process noise w and measurement noise v are Gaussian white
%   noises with covariance:
%
%         E ([w;v] * [w',v']) = QWV.
%
%   LQG uses LQR and KALMAN to compute the LQG regulator. The state-space
%   model SYS should specify the A,B,C,D matrices, see SS for details.
%   LQG can be used for both continuous- and discrete-time plants. By
%   default, LQG uses x[n|n-1] as state estimate in discrete time, see
%   KALMAN for details. To use x[n|n] instead, type
%      KLQG = LQG(SYS,QXU,QWV,...,'current')
%
%   KLQG = LQG(SYS,QXU,QWV,QI) computes an LQG servo-controller KLQG that
%   uses the setpoint command r and measurements y to generate the control
%   signal u. KLQG has integral action to ensure that the output y tracks
%   the command r.
%
%                                       | w             | v
%                       .----------.    |   .-------.   |
%          r   -------->|          |    '-->|       |   V
%                       |   KLQG   |------->|  SYS  |---O--.-----> y
%          y   .------->|          |  u     |       |      |
%              |        '----------'        '-------'      |
%              |                                           |
%              '-------------------------------------------'
%
%   The LQG servo-controller minimizes the cost function
%
%         J(u) = Integral [x',u'] * QXU * [x;u] + xi' * Qi * xi dt
%
%   where xi is the integral of the tracking error r-y and SYS,w,v are as
%   described above. For MIMO systems, r, y, and xi must have the same
%   length. LQG uses the commands LQI and KALMAN to compute KLQG.
%
%   KLQG = LQG(SYS,QXU,QWV,QI,'1dof') computes a one-degree-of-freedom
%   servo-controller that takes e=r-y rather than [r;y] as input.
%
%   KLQG = LQG(SYS,QXU,QWV,QI,'2dof') is equivalent to LQG(SYS,QXU,QWV,QI)
%   and produces the two-degree-of-freedom servo-controller shown above.
%
%   [KLQG,INFO] = LQG(SYS,QXU,QWV,...) returns the controller and estimator
%   gains in the structure INFO. The controller equations are as follows:
%
%   Continuous:
%      dx_e = A x_e + B u + L (y - C x_e - D u)
%         u = - Kx x_e - Ki xi
%
%   Discrete:
%      x[n+1|n] = A x[n|n-1] + B u[n] + L (y[n] - C x[n|n-1] - D u[n])
%      delayed:  u[n] = - Kx x[n|n-1] - Ki xi[n]
%      current:  u[n] = - Kx x[n|n] - Ki xi[n] - Kw w[n|n]
%                     = - Kx x[n|n-1] - Ki xi[n] - (Kx*Mx+Kw*Mw) innov[n]
%                with innov[n] = y[n] - C x[n|n-1] - D u[n].
%
%   See also TVLQR, CDRE, TVKALMAN

% Extract System Matrices
[A,B,C,D] = ssdata(SYS); %#ok<ASGLU>

% Default Horizon
T0 = A.Time(1);
Tf = A.Time(end);
Tspan = [T0,Tf];

% Design LQR controller over a finite horizon
Kc = tvlqr(A,B,Q,R,[],[],Tspan);

% Design Kalman Filter over a finite horizon
Kf = tvkalman(SYS,Qn,Rn,Fn,Tspan);

% Augment the results to build final controller
Klqg = tvss(A-B*Kc-Kf*C,L,-Kc,0);



