function out = tvsens(G,yBar,pBar,Rho,Nt,Np)
%% SENS Performs first order sensitivity analysis
%
% Inputs:
%
% G     := LTV system (Specified as TVSS)                   [Required]
% yBar  := Trim outputs of interest (Specified as TVMAT)    [Required]
% pBar  := Nominal Parameter Array (Np x 1)                 [Required]
% Rho   := Percentage variation for parameters (Np x 2)     [Required]
% Nt    := Number of time grid points for analysis          [Optional]
% Np    := Number of parameters to consider for analysis 	[Optional]
%
% Outputs:
%
% out   := Structure include bounds on trajectories at analysis grid points

%% Input Processing
% Check number of input and output arguments
narginchk(4,6);
nin = nargin;
[Ny,Nu] = size(G);
switch nin
    case 4
        % Default 20 points equaly spaced on the horizon
        Nt = 20;
        % Consider all the parameters that are in the input
        Np = Nu;
    case 5
        Np = Nu;
end

% Check nominal parameter vector pBar and percentage variation Rho
lpBar = length(pBar);
[szr1,szr2] = size(Rho);
if (szr1~=lpBar) || (szr2~=2)
    error('Parameter array must be Nx1 and percentage variation about nominal must be specified through Nx2 array.')
end

% Time grid for analysis
[T0,Tf] = getHorizon(G);
Tgrid = linspace(T0,Tf,Nt);

% XXX: Check that Ybar is of dimention 1
[mYbar,nYbar] = size(yBar);
if mYbar~=1 || nYbar~=1
    error('Size of Ybar must be 1x1.')
end

%% Simulate Linear System
% States are perturbation around the nominal trajectories

% Propogate Individual Sensitivities
U = evalt(tvmat(eye(Nu)),G.Time);
delY = cell(Np,1);
for i = 1:Np
    delY{i} = tvlsim(G,U(:,i));
end

%% Linear programming
% Solve optimization problem with box constraints at each time instant
YUpp = zeros(Ny,1,Nt);
YLow = zeros(Ny,1,Nt);
PLow = zeros(Np,1,Nt);
PUpp = zeros(Np,1,Nt);
tComp = zeros(Nt,1);

% Optimization Options
opts = optimoptions('linprog','Display','off');

% Outer for loop
for i = 2:Nt
    %% Output at time t on the Tgrid
    Yt = tvsubs(yBar,Tgrid(i));
    delYt = cellfun(@(x) tvsubs(x,Tgrid(i)),delY);
    
    %% Linear Programming
    
    % Create Optimization Variable
    % Assumes that Rho(:,1) is lower bound and Rho(:,2) is upper bound
    ub = zeros(Np,1);
    lb = zeros(Np,1);
    for j = 1:Np
        lb(j) = sign(pBar(j))*(pBar(j)*Rho(j,1)/100);
        ub(j) = sign(pBar(j))*(pBar(j)*Rho(j,2)/100);
    end
    du = optimvar('du',Np,'LowerBound',lb,'UpperBound',ub);
    
    % Lower and Upper Bound Problems
    LBProb = optimproblem('Objective',Yt+delYt'*du,'ObjectiveSense','min');
    UBProb = optimproblem('Objective',Yt+delYt'*du,'ObjectiveSense','max');
    
    % Solve LP
    t0 = tic;
    [sol1,YLow(:,:,i),ef1,op1] = solve(LBProb,'options',opts); %#ok<ASGLU>
    [sol2,YUpp(:,:,i),ef2,op2] = solve(UBProb,'options',opts); %#ok<ASGLU>
    tComp(i) = toc(t0);
    
    PLow(:,:,i) = sol1.du;
    PUpp(:,:,i) = sol2.du;
end

%% Output structure
out = struct;
out.PLow = tvmat(PLow,Tgrid);
out.PUpp = tvmat(PUpp,Tgrid);
out.YLow = tvmat(YLow,Tgrid);
out.YUpp = tvmat(YUpp,Tgrid);
out.tComp = tComp(2:Nt);
end