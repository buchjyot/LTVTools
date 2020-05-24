function cost = QSRFCost(X,U,Q,R,varargin)
%% QSRFCost
% This function computes the finite horizon integral qudratic cost using
% the specified cost matrices Q, R, S, and F.
%
% cost = X(T)'*F*X(T) + int_0^T [X(t);U(t)]' [Q S;S' R] [X(t);U(t)] dt

%% Input Processing
narginchk(4,6);
nin = nargin;

% Get Horizon
[~,Tf] = getHorizon(X);

% Cost Matrices
S = [];
F = [];
switch nin
    case 5
        S = varargin{1};
    case 6
        S = varargin{1};
        F = varargin{2};
        XTf = tvsubs(X,Tf);
end

%% Compute Cost
integrand = [X;U]'*[Q S;S' R]*[X;U];
integralCost = trapz(integrand.Time,integrand.Data);
if nin >= 6
    terminalCost = XTf'*F*XTf;
    cost = terminalCost + integralCost;
else
    cost = integralCost;
end
