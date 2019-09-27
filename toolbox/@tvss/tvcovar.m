function [P,Q,info] = tvcovar(SYS,varargin)
% TVCOVAR calculates the stationary covariance of the output y of finite
% horizon LTV model sys driven by Gaussian white noise inputs w.
%
% P = tvcovar(sys,W) returns the steady-state output response covariance
% P = E(y*y') given the (optional) noise intensity W,
%
%     E(w(tau)*w(tau)') = W*delta(t-tau)
%
% NOTE: By default W = 1
%
% P = tvcovar(sys,W,OPT) returns the steady-state output response covariance
% P = E(y*y') given the noise intensity W, uses the OPT as tvOdeOptions
%
% Following lyapunov diffrential equation (CDLE) is solved
% Qdot(t) = A(t)*Q(t) + Q(t)*A'(t) + B(t)*W(t)*B'(t), Q(0) = 0
%
% Assumes 0 initial covariance

%% Input Processing
narginchk(1,3);
nin = nargin;

% Obtain default Opt
switch nin
    case 1
        Opt = tvodeOptions;
        W = 1;
    case 2
        if isa(varargin{1},'tvodeOptions')
            Opt = varargin{1};
            W = 1;
        else
            Opt = tvodeOptions;
            W = varargin{1};
        end
    case 3
        W = varargin{1};
        Opt = varargin{2};
end

% Check Feedthrough condition
[A,B,C,D] = ssdata(SYS);
if any(D.Data(:))
    warning('ltvtools:tvcovar:nonzerofeedthrough','The feedthrough must be zero.');
    P = [];
    Q = [];
    return;
end

% Check if horizon is finite
ltvutil.verifyFH(SYS);

% Compute Controllability Gramian of modified system
SYSm = tvss(A,B*sqrt(W),C,D);
[Q,~,info]  = tvgram(SYSm,'c',Opt);
[C,Q] = evalt(C,Q,union(C.Time,Q.Time));
P = C*Q*C';
end