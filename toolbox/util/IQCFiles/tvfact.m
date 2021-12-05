function G = tvfact(varargin)
% TVFACT  Spectral factorization of linear systems on finite horizons.

%   G = SPECTRALFACT(F,R) handles the case when H is specified in
%   factored form as H = F' * R * F. The spectral factorization is then
%   computed without explicitly forming H.

%% Input Checking
narginchk(3,4)
nin = nargin;
switch nin
    case 3
        Tspan = varargin{3};
        Opt = tvodeOptions;
    case 4
        Tspan = varargin{3};
        Opt = varargin{4};
end

% System Data
sys = varargin{1};
[A,B,C,D] = ssdata(sys);
[Nx,~] = size(B);

%% Build Cost Matrices
M11 = varargin{2};
Qc = C'*M11*C;
Sc = C'*M11*D;
Rc = D'*M11*D;
Fc = zeros(Nx);
W = chol(Rc);

%% Solve Riccati Equation
Tspan = flip(sort(Tspan));
X = cdre(tvmat(A),B,Qc,Rc,Sc,[],Fc,Tspan,Opt);

%% Final Output
[A,B,Sc] = evalt(tvmat(A),B,Sc,X.Time);
G = tvss(A,B,W\(X*B + Sc)',W);
end