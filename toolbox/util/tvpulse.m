function out = tvpulse(T0,Tf,t1,t2,M1,M2,varargin)
%% TVPULSE Creates a pulse of desired magnitude
%
% Inputs:
%
% [T0,Tf] : Horizon
% [t1,t2] : Pulse Duration
% [M1,M2] : Transition Level
%  Ts     : Time grid parameter
%
% Output:
%
% out : tvmat that is defined on [T0,Tf], and have a transition from level
% M1 to M2 such that it has magnitude M2 during [t1,t2] time duration.
%
%       out = M1, T0<=t<=t1
%           = M2, t1<=t<=t2
%           = M1, t2<=t<=Tf
%
% Syntax:
%
% >> A = tvpulse(0,10,2,4,0,0.1,0.1)

narginchk(6,7);
nin = nargin;

% Default
Ts = 0.01;

switch nin
    case 7
        Ts = varargin{1};
end

warning('off','ltvtools:evalt:extrapolate');
A = tvmat(M2);
At = evalt(A,[t1,t2]);
out = evalt(At,T0:Ts:Tf);
out.Data(out.Data == 0) = M1;
warning('on','ltvtools:evalt:extrapolate');
end