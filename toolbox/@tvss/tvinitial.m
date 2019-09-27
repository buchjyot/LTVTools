function [Y,X] = tvinitial(G,x0,Opt)
% TVINITIAL Simulate the Initial Condition Response of LTV System G
%
% [Y,X] = tvinitial(G,X0) where G is a TVSS and X0 is initial condition

% Input Processing
narginchk(1,3);
nin = nargin;

% Get sizes of states and inputs
[A,B] = ssdata(G);
Nx = size(A,1);
Nu = size(B,2);

% Unforced Response
U = tvmat(zeros(Nu,1,size(G.Time,1)),G.Time);

switch nin
    case 1
        % Assume Zero Initial Conditions if nothing is specified
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
        warning('Initial conditions were not specified, assuming zero initial conditions.');
    case 2
        Opt = tvodeOptions;
end

[Y,X] = tvlsim(G,U,x0,Opt);