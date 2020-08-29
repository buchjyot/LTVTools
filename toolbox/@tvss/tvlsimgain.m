function [g,Y,X] = tvlsimgain(G,U,NE,x0,Opt)
%% TVSIMGAIN
%
% User specify an LTV system G and a candidate bad input U. We simulate the
% system G using tvlsim and compute the induced gain for nominal
% performance specified using NE.
%
% Optional arguments:
% NE is default to 0 (L2toL2) and Opt is default to tvodeOptions
%
% NOTE:
% Let G be a TVSS system with NY by NU IO dimentions
% NE is an optional argument with default set to 0.
% NE = NY means 'L2toE' Norm
% 0 < NE < NY then we have 'L2toE' and 'L2toL2' mixed cost terms.

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
Nx = order(G);
NY = size(G,1);

% Input Processing
narginchk(2,5);
nin = nargin;
switch nin
    case 2
        NE = 0;
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
    case 3
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
    case 4
        Opt = tvodeOptions;
end

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY.');
end

% Linear simulation
[Y,X] = tvlsim(G,U,[T0,Tf],x0,Opt);

% Compute induced gain
NL2 = NY-NE;
YL2 = Y(1:NL2,1);
YE  = tvsubs(Y(NL2+1:end,1),Tf);
g   = (sqrt(norm(YE)^2 + tvnorm(YL2)^2))/tvnorm(U);
end