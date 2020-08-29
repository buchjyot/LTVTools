function [Y,X] = tvimpulse(G,Tf,Opt)
%% TVIMPULSE
% Computes the impulse response of the LTV system on finite time horizon. 
%
% Example:
%  
%   >> Gss = ss(-1,1,1,0);
%   >> G = tvss(Gss,0:0.1:10);
%   >> [y,x] = tvimpulse(G,5)

% Input Processing
narginchk(1,3);
nin = nargin;
GTime = G.Time;
[GT0,GTf] = getHorizon(G);

T0 = GT0;
switch nin
    case 1
        Tf = GTf;
        Opt = tvodeOptions;
    case 2
        Opt = tvodeOptions;
end

% Assume Zero Input
[Ny,Nu] = size(G);
U = tvmat(zeros(1,1,size(GTime,1)),GTime);

% Memory Allocation
Y = cell(Ny,Nu);
X = cell(Ny,Nu);

% Use the fact that impulse response is an initial condition response with
% x0 set to the B(0) vector.
for i = 1:Ny
    for j = 1:Nu
        Gi = tvss(G.Data(i,j,:),GTime,G.InterpolationMethod);
        [~,Bi] = ssdata(Gi);
        x0 = tvsubs(Bi,T0);
        [Y{i,j},X{i,j}] = tvlsim(Gi,U,[T0,Tf],x0,Opt);
    end
end
end