function [Y,X] = tvinitial(G,x0,Tf,Opt)
%% TVINITIAL Simulate the Initial Condition Response of LTV System G
%
% [Y,X] = tvinitial(G,X0) where G is a TVSS and X0 is initial condition
%
% [Y,X] = tvinitial(G,X0,Tf,Opt) where Tf is a final time and Opt is
% tvodeOptions

% Input Processing
narginchk(1,3);
nin = nargin;

% Get sizes of states and inputs
Nx = order(G);
[~,Nu] = size(G);

% Unforced Response
GTime = G.Time;
GT0 = GTime(1);
GTf = GTime(end);
U = tvmat(zeros(Nu,1,size(GTime,1)),GTime);
T0 = GT0;

switch nin
    case 1
        Tf = GTf;
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
        warning('Initial conditions were not specified, assuming zero initial conditions.');
    case 2
        Tf = GTf;
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
        warning('Initial conditions were not specified, assuming zero initial conditions.');
    case 3
        Opt = tvodeOptions;
end

% Call TVLSIM
[Y,X] = tvlsim(G,U,[T0,Tf],x0,Opt);