function [wcg,wcinfo] = tvwcgain(G,Delta,varargin)
%% TVWCGAIN Worst-Case Gain of uncertain LTV system on a finite horizon
%
% Inputs
%   G - Nominal LTV system
%   Delta - Uncertainty charecterized by udyn
%   NE - Number of outputs penalized in Euclidean sense
%   Opt - tvwcOptions set for tvwcgain
%
% Outputs
%   wcg - Gain Upper bound
%   info - Structure of solution info
%
% NOTE: RDE/DLMI based analysis will be performed on the horizon, in which
% the LTV system is defined.

%% Input Processing
nin = nargin;
narginchk(2,4);

NE = []; Opt = [];
switch nin
    case 3
        if isa(varargin{1},'tvwcOptions')
            Opt = varargin{1};
        elseif isa(varargin{1},'double')
            NE = varargin{1};
        end
    case 4
        NE = varargin{1};
        Opt = varargin{2};
end

% Use Default Options
if isempty(Opt)
    Opt = tvwcOptions;
end

% Time Horizon
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

%% Read Options
DispFlag = isequal(Opt.Display,'on');
StopTol = Opt.StopTol;
LMIopt = Opt.LMIOptions;
RDEopt = Opt.RDEOptions;

% Maximum # of iterations
Niter = Opt.MaxIter;

% Time Grid for LMI Approach
tlmi = linspace(T0,Tf,Opt.Nlmi);
tlmi = tlmi(:);

% Time Points for Spline Basis Functions
Nsp = Opt.Nsp;
tSp = linspace(T0,Tf,Nsp);

% Create Spline Basis Functions
Ps = tvmat(reshape(eye(Nsp),Nsp,1,Nsp),tSp,'Spline');

%% Memory Allocation
glmi = zeros(Niter,1);
grde = zeros(Niter,1);
AllIter(Niter,1) = struct('TotalTime',[],'LMI',[],'RDE',[]);

%% Read UserData
[Nw,Nv] = size(Delta);
[NY,NU] = size(G);

%% Plant Scalling Based Specified ULevel
% Check uncertainty level and scale the plant accordingly
sqrtUL = sqrt(Opt.ULevel);
if ~isequal(sqrtUL,1)
    G = blkdiag(sqrtUL,eye(NY-Nv))*G*blkdiag(sqrtUL,eye(NU-Nw));
end

%% Iterations

% Display starting
if DispFlag
    fprintf(' Starting IQC upper bound iterations...\n');
end
    
% Begin Timing
t0 = tic;

for i=1:Niter
    
    if DispFlag
        fprintf(' Iter# = %d: ',i);
    end
    
    %% Evaluate data on LMI time grid
    G1 = evalt(G,tlmi);
    Ps1 = evalt(Ps,tlmi);
    Psdot1 = tvdiff(Ps,tlmi);
    if i>1 && isfinite(grde(i-1))
        % Use RDE solution from previous iteration (if grde<inf)
        [Pm,Pmdot] = evalRDE(solrdePrev,tlmi);
    else
        Pm = tvmat; Pmdot = Pm;
    end
    
    %% Finite Horizon: LMI Condition
    [glmi(i),LMIinfo] = ltvutil.fhlmiEngine(G1,Delta,Ps1,Psdot1,Pm,Pmdot,LMIopt,NE);
    if DispFlag
        fprintf(' DLMI Gain Bound = %4.4f,',glmi(i));
    end
    
    % Evaluate LMI solution (for debugging) to show eLMI1<=0
    % [LMI1,eLMI1] = evalLMI(G1,Delta,LMIinfo.P,LMIinfo.Pdot,glmi(i),glmi(i),LMIinfo,NE);
    % figure; tvplot(eLMI1);
    
    %% Finite Horizon: RDE + Bisection
    [gbnds,RDEinfo] = ltvutil.fhrdeEngine(G,Delta,LMIinfo,RDEopt,NE);
    
    % Choose an Upper Bound
    grde(i) = gbnds(2);
    if isfinite(grde(i))
        trde = fliplr(RDEinfo.Upper.sol.x);
        solrde = RDEinfo.Upper.sol;
    else
        trde = fliplr(RDEinfo.Lower.sol.x);
        solrde = RDEinfo.Lower.sol;
    end
    
    if DispFlag
        fprintf(' RDE Gain Bound = %4.4f,',grde(i));
        if isinf(grde(i))
            fprintf(' trde = %4.3f',trde(1));
        end
    end
    
    %% Termination Condition: LMI and RDE gains are "close"
    if abs(grde(i) - glmi(i)) < StopTol*glmi(i)
        break;
    end
    
    %% Compare LMI and RDE solutions and update time grid
    if glmi(i) < grde(i)
        % LMI cost was better indicating the LMI time grid was too coarse
        
        % Create dense time grid
        tDense = union(trde,tlmi);
        t1 = trde(1);
        if t1>0
            % RDE solution diverged at t1. Add time points in [0,t1]
            tdensity = numel( tDense(tDense>=t1) ) / (Tf-t1);
            tadd = linspace(0,t1, ceil(t1*tdensity) );
            tDense = union(tDense,tadd);
        end
        
        % Evaluate LMI using on dense time grid
        G2 = evalt(G,tDense);
        Ps2 = evalt(Ps,tDense);
        Psdot2 = tvdiff(Ps,tDense);
        if i>1 && isfinite(grde(i-1))
            % Convergent RDE solution on previous iteration
            [Pm2,Pmdot2] = evalRDE(solrdePrev,tDense);
        else
            % No convergent RDE solution on previous iteration
            Pm2 = tvmat; Pmdot2 = Pm2;
        end
        [Pval,Pdotval] = evalP(Ps2,Psdot2,Pm2,Pmdot2,LMIinfo.Y,LMIinfo.lamP);
        [~,eLMI2] = evalLMI(G2,Delta,Pval,Pdotval,glmi(i),LMIinfo,NE);
        
        % Update the LMI time grid.
        eLMI2 = eLMI2.Data(:);
        emax = max(eLMI2);
        addidx = [];
        for j=1:numel(tlmi)-1
            idx = find( tDense>tlmi(j) & tDense<tlmi(j+1) );
            if ~isempty(idx)
                [~,idx2] = max( eLMI2(idx) );
                if eLMI2(idx(idx2))> 0.25*emax
                    % We do not know before hand that howmany grid points
                    % will be added so can not pre-allocate memory.
                    addidx = [addidx; idx(idx2)]; %#ok<AGROW>
                end
            end
        end
        tadd = tDense(addidx);
        
        if DispFlag
            fprintf(' Adding t = ');
            for j=1:numel(tadd)
                fprintf('%4.3f, ',tadd(j));
            end
        end
        
        % figure(3)
        % plot(tlmi,eLMI1,'b',tlmi,eLMI1,'bx',tt,eLMI2,'r',...
        %       tadd,eLMI2(addidx),'rx')
        % ylim(max(eLMI2(:))*[-1 1]);
        % drawnow;
        
        tlmiNew = sort([tlmi; tadd]);
        if isequal(tlmiNew,tlmi)
            % Terminating - Same Time Grid
            break;
        end
        tlmi = tlmiNew;
    end
    
    % Store Iteration Info
    solrdePrev = solrde;
    AllIter(i) = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
    
    % Newline for next iteration
    if DispFlag && i~=Niter
        fprintf(newline);
    end
end

%% Process Outputs
% Store Iteration Info (if iteration terminated before Niter)
if i < Niter
    AllIter(i) = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
end

% Display logs
if DispFlag
    if i == Niter
        fprintf('\n Maximum number of iteration reached.');
    else
        fprintf('\n Stopping tolerance satisfied: terminating analysis iterations.');
    end
end

% Store outputs
[wcg,idx] = min(grde(1:i)); % grde(i);
AllIter = AllIter(1:i);
tTotal = toc(t0);
if DispFlag
    fprintf('\n Final Results:');
    fprintf(' RobustPerfUB = %4.4f,',wcg);
    fprintf(' TotalCompTime = %4.4f\n',tTotal);
end
wcinfo = AllIter(idx);
wcinfo.AllIter   = AllIter;
wcinfo.TotalIter = i; 
wcinfo.TotalTime = tTotal;
% wcinfo.WCIter    = idx;
end

%% LOCAL Function: evalP
function [P,Pdot] = evalP(Ps,Psdot,Pm,Pmdot,Y,lamP)
% Combine scalar and matrix bases functions
%    P = ps1*Y1+ ... psN*YN + lamP*Pm
% The derivative Pdot is computed simiarly.

% Sum terms with scalar basis functions
Ps = Ps(:);
Ns = size(Ps,1);
P = Ps(1)*Y(:,:,1);
Pdot = Psdot(1)*Y(:,:,1);
for j=2:Ns
    P = P + Ps(j)*Y(:,:,j);
    Pdot = Pdot + Psdot(j)*Y(:,:,j);
end

% Include matrix basis function (if it exists)
if ~isempty(Pm)
    % XXX
    Pm.InterpolationMethod = 'Spline';
    Pmdot.InterpolationMethod = 'Spline';
    
    P = P + lamP*Pm;
    Pdot = Pdot + lamP*Pmdot;
end
end

%% LOCAL Function: evalRDE
function [P,Pdot] = evalRDE(solrde,t)
% Evaluate RDE solution on a specified time grid.

% Dimensions
Nt = numel(t);
Nx = round(sqrt(  size(solrde.y,1) ));

% Using deval to evaluate Pdot is different (less accurate?) than
% directly calling the ODEFH to compute the derivative (as done below)
%     [P,Pdot] = deval(solrde,t);
%     P = reshape( P, [Nx Nx Nt]);
%     Pdot = reshape( Pdot, [Nx Nx Nt]);

% Construct P
P = deval(solrde,t);
P = reshape(P,[Nx Nx Nt]);

% Evaluate Pdot
odefh = solrde.extdata.odefun;
Pdot = zeros(Nx,Nx,Nt);
for i=1:Nt
    Pdot(:,:,i)=reshape( odefh(t(i),P(:,:,i)), [Nx Nx]);
end

% Ensure symmetry
for i=1:Nt
    P(:,:,i) = ( P(:,:,i)+P(:,:,i)' )/2;
    Pdot(:,:,i) = ( Pdot(:,:,i)+Pdot(:,:,i)' )/2;
end

% Return as TVMAT
P = tvmat(P,t);
Pdot = tvmat(Pdot,t);
end

%% LOCAL Function: evalLMI
function [LMI,eLMI] = evalLMI(G,Delta,P,Pdot,g,LMIinfo,NE)
% This function evaluates the finite-horizon gain condition on a time grid
% for given values of the IQC, storage function, and gain. **The function
% assumes G is LTV and Delta is a SISO, unit norm bounded, LTI uncertainty.
%
% NOTE - All data (G,P,Pdot) must have the same time grid t

% Only required for time-varying IQC multiplier lambda
[Nw,Nv] = size(Delta);
[Type,Psi] = ltvutil.iqcEngine(Delta);
tDense = G.Time;
if isequal(Type,1)
    X11 = evalt(LMIinfo.X11,tDense); %#ok<*UNRCH>
else
    X11 = LMIinfo.X11;
end
X = blkdiag(X11,-X11);

% Form Extended System
[AAe,BBe,CCe,DDe,~,Nz,~,Nd] = ltvutil.extsystem(G,Psi,Nv,Nw,NE);
n = size(AAe,1);
CC1e = CCe(1:Nz,:);
DD1e = DDe(1:Nz,:);
CC2e = CCe(Nz+1:end,:); % Should be zero for pure Euclidean Cost
DD2e = DDe(Nz+1:end,:); % Should be zero for pure Euclidean Cost

% Evaluate LMI
P.InterpolationMethod = 'Linear';
Pdot.InterpolationMethod = 'Linear';
gsq = g^2;
LMI = [AAe'*P+P*AAe+Pdot P*BBe; BBe'*P zeros(Nd+1)] ...
    + blkdiag( zeros(n+1), -gsq*eye(Nd) ) + [CC1e DD1e]'*X*[CC1e DD1e] ...
    + [CC2e DD2e]'*[CC2e DD2e];

% Evaluate max real eigenvalue
eLMI = max(real( eig(LMI) ));
end