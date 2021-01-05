function [bnd,dwc,info] = tvnormaff(sysAug,varargin)
%% tvnormaff
% Given system G, uncertainty Delta, trim input vbar, and norm of the
% disturbance alpha, this function computes the Euclidean norm bound on the
% output at the final time, once Delta is wrapped in to the system G.
%
% Optional arguments
% alpha is default to 1
% Opt is specified using tvnormOptions
%
% Syntax:
% [bnd,dwc,info] = tvnormaff(sysAug,alpha,NE,Opt);

% Minimum 3 intput arguments are required, maximum 4.
narginchk(1,4);
nin = nargin;
nout = nargout;

% Default
alpha = [];
NE = [];
Opt = tvnormOptions;
switch nin
    case 2
        alpha = varargin{1};
        
    case 3
        alpha = varargin{1};
        NE = varargin{2};
        
    case 4
        alpha = varargin{1};
        NE = varargin{2};
        Opt = varargin{3};
        
end
if isempty(NE)
    NE = size(sysAug,1);
end
if isempty(alpha)
    alpha = 1;
end
cdleOpt = tvodeOptions('OdeSolver',Opt.OdeSolver,'OdeOptions',Opt.OdeOptions);
Display = isequal(Opt.Display,'on');

% Verify G is finite horizon
ltvutil.verifyFH(sysAug);
[T0,Tf] = getHorizon(sysAug);

% 09/25/2020 Make sure NY == NE
% i.e. we do not support induced L2 penalties on output at this moment
NY = size(sysAug,1);
if ~isequal(NY,NE)
    error('NY must be equal to NE.');
end

% State-space data
[Abar,Bbar,Cbar,Dbar] = ssdata(sysAug);
NL2 = NY-NE;
sysAdj = tvss(-Abar',-Cbar(1:NL2,:)',Bbar',Dbar(1:NL2,:)');

% Make sure the feedthrough is 0.
DTf = tvsubs(Dbar,Tf);
if any(DTf(NL2+1:end,:))
   error('Feedthrough term d->e must be zero for well-posed Euclidean penalty'); 
end

%% Compute Output Controllability Gramian
CTbar           = tvsubs(Cbar,Tf);
[Pbar,~,Pinfo]  = cdle(Abar,Bbar,zeros(order(sysAug)),[T0 Tf],cdleOpt); % Pbar  = tvgram(sysAug,'c',Opt);
PTbar           = tvsubs(Pbar,Tf);
WTbar           = CTbar*PTbar*CTbar';

% Eigenvalue decomposition
[Qtmp,Dtmp] = eig(WTbar);

% Sort eigenvalues
[lambdaVec,ind] = sort(diag(Dtmp),'descend');
Q = Qtmp(:,ind);

%% Compute gammaTilde
% simulate unforced system
x0 = sysAug.UserData;
if isempty(x0)
    warning(' UserData is empty, simulating system with zero initial condition.');
end
[Y,~] = tvlsim(sysAug,[],[T0,Tf],x0);

% compute gamma
gamma = tvsubs(Y,Tf);
gammaTilde = Q'*gamma;

%% Algorithm
ne = numel(gammaTilde);
lambda1 = lambdaVec(1);

% Compute alphaMax
epsilon = 1e-6;
rtmp = find(lambdaVec > lambda1 - epsilon);
r = rtmp(end);

if all(gammaTilde(1:r) == 0)
    if Display
        fprintf(' alphaMax finite, ');
    end
    alphaMaxSq = sum(gammaTilde(r+1:end).^2.*lambdaVec(r+1:end)./(lambda1-lambdaVec(r+1:end)).^2);
else
    if Display
        fprintf(' alphaMax infinite, ');
    end
    alphaMaxSq = Inf;
end

% Compute zStar, eMaxSq
if alpha^2 >= alphaMaxSq
    if Display
        fprintf('Special Case\n');
    end
    
    % compute beta
    beta = sqrt((alpha^2-alphaMaxSq)/lambda1);
    
    % compute zStarR
    LambdaR = diag(lambdaVec(r+1:end));
    gammaTildeR = gammaTilde(r+1:end);
    zStarR = (lambda1*eye(ne-r)-LambdaR)\gammaTildeR;
    
    % compute zStar
    e1 = zeros(r,1); e1(1) = 1;
    zStar = [beta*e1;zStarR];
    
    % compute eMaxSq
    eMaxSq = lambda1^2*beta^2 + lambda1^2*sum(gammaTilde(r+1:end).^2./(lambda1-lambdaVec(r+1:end)).^2);
    
else
    if Display
        fprintf('Regular Case\n');
    end
    
    % compute muAlpha
    f = @(mu) alpha^2 - sum(gammaTilde.^2.*lambdaVec./(mu-lambdaVec).^2);
    muRng = [1.0000001, 1e5]*lambda1;
    muAlpha = fzero(f,muRng);
    
    % compute zStar
    zStar = (muAlpha*eye(ne)-diag(lambdaVec))\gammaTilde;
    
    % compute eMaxSq
    eMaxSq = muAlpha^2*sum(gammaTilde.^2./(muAlpha-lambdaVec).^2);
end

bnd = sqrt(eMaxSq);

%% Construct Disturbance
if nout > 1
    adjBCNonzero = CTbar'*Q*zStar;
    [dwc, ~] = tvlsim(sysAdj,[],[Tf,T0],adjBCNonzero);
end
info = [];
info.Pinfo = Pinfo;
end