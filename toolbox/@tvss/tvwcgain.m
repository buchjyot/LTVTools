function [wcg,info] = tvwcgain(G,varargin)
%% TVWCGAIN Worst-Case Gain of uncertain LTV system on a finite horizon
%
% Inputs
%   G - Nominal LTV system (as a TVSS). The uncertain system is
%       Fu(G,Delta) where Delta is a SISO, uncertainty
%       satisyfing the unit norm bound ||Delta|| <= 1.
%   NE - Number of outputs penalized in Euclidean sense
%   tvwcOptions - Options set for tvwcgain
%
% Outputs
%   wcg - Gain Upper bound
%   info - Structure of solution info
%
% NOTE: RDE/DLMI based analysis will be performed on the horizon, in which
% the LTV system is defined.

%% Input Processing
nin = nargin;
% nout = nargout;
narginchk(1,3);

NE = []; Opt = [];
switch nin
    case 2
        if isa(varargin{1},'tvwcOptions')
            Opt = varargin{1};
        elseif isa(varargin{1},'double')
            NE = varargin{1};
        end
    case 3
        NE = varargin{1};
        Opt = varargin{2};
end

% Use Default Options
if isempty(Opt)
    Opt = tvwcOptions;
end

% Time Horizon
ltvutil.verifyFH(G);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% XXX Check that Euclidean penalty is valid

%% Read Options
DispFlag = isequal(Opt.Display,'on');
StopTol = Opt.StopTol;
LMIopt = Opt.LMIOptions;
RDEopt = Opt.RDEOptions;

% Maximum # of iterations
Niter = Opt.MaxIter;

% Time Horizon
[T0,Tf] = getHorizon(G);

% Time Grid for LMI Approach
tlmi = linspace(T0,Tf,Opt.Nlmi);
tlmi = tlmi(:);

% Time Points for Spline Basis Functions
Nsp = Opt.Nsp;
tSp = linspace(T0,Tf,Nsp);

% Create Spline Basis Functions
Ps = tvmat(reshape(eye(Nsp),Nsp,1,Nsp),tSp,'Spline');

% See user data for IQC Parameterization
IQCParam = G.UserData;
if isempty(IQCParam)
    warning('ltvtools:tvwcgain:usingDefaultIQCParam','UserData property is empty, assuming default IQC parameterization v = p = 0');
    v = 0;p = 0;
else
    v = IQCParam.v;
    p = IQCParam.p;
end

%% Memory Allocation
glmi = zeros(Niter,1);
grde = zeros(Niter,1);
info = cell(Niter,1);

%% Plant Scalling Based Specified ULevel
% XXX: This will be part of tvuss later.
Nv = 1;
Nw = 1;

% Check uncertainty level and scale the plant accordingly
sqrtUL = sqrt(Opt.ULevel);
[NY,NU] = size(G);
if ~isequal(sqrtUL,1)
    G = blkdiag(sqrtUL,eye(NY-Nv))*G*blkdiag(sqrtUL,eye(NU-Nw));
end

% Flag whether you want to consider timevarying IQC matrix
tvIQCFlag = Opt.TimeVaryingIQC;

% Choose LMI Function based on tvIQCFlag
if tvIQCFlag
    FHLMIFCN = @ltvutil.tvfhlmi;
else
    FHLMIFCN = @ltvutil.fhlmi;
end

%% Iterations
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
    [glmi(i),X11,Plmi,Pdotlmi,Y,lamP] = ...
        FHLMIFCN(G1,v,p,Ps1,Psdot1,Pm,Pmdot,LMIopt,NE);
    if DispFlag
        fprintf('LMI Gain Bound = %4.4f,',glmi(i));
    end
    LMIinfo = struct('X11',X11,'P',Plmi,'Pdot',Pdotlmi,'Y',Y,'lamP',lamP);
    
    % Evaluate LMI solution (for debugging) to show eLMI1<=0
    % [LMI1,eLMI1] = evalLMI(G1,NE,v,p,glmi(i),X11,Plmi,Pdotlmi);
    % figure; tvplot(eLMI1);
    
    %% Finite Horizon: RDE + Bisection
    if tvIQCFlag
        X11RDE = evalt(tvmat(X11),G.Time);
    else
        X11RDE = X11;
    end
    [gbnds,RDEinfo] = ltvutil.fhrde(G,v,p,Tf,X11RDE,[],RDEopt,NE);
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
        
        % Evaluate LMI using (P,Pdot,X11,g) but on dense time grid
        G2 = evalt(G,tDense);
        if tvIQCFlag
            X11Dense = evalt(X11,tDense); %#ok<*UNRCH>
        else
            X11Dense = X11;
        end
        Ps2 = evalt(Ps,tDense);
        Psdot2 = tvdiff(Ps,tDense);
        if i>1 && isfinite(grde(i-1))
            % Convergent RDE solution on previous iteration
            [Pm2,Pmdot2] = evalRDE(solrdePrev,tDense);
        else
            % No convergent RDE solution on previous iteration
            Pm2 = tvmat; Pmdot2 = Pm2;
        end
        
        [Pval,Pdotval] = evalP(Ps2,Psdot2,Pm2,Pmdot2,Y,lamP);
        [~,eLMI2] = evalLMI(G2,NE,v,p,glmi(i),X11Dense,Pval,Pdotval);
        
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
    info{i} = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
    
    % Newline for next iteration
    if DispFlag && i~=Niter
        fprintf(newline);
    end
end

%% Process Outputs
% Store Iteration Info (if iteration terminated before Niter)
if i < Niter
    info{i} = struct('TotalTime',toc(t0),'LMI',LMIinfo,'RDE',RDEinfo);
end

% Store outputs
[wcg] = min(grde(1:i)); % grde(i);
info = info(1:i);
tTotal = toc(t0);
if DispFlag
    fprintf('\n Final Results:');
    fprintf(' RobustGain = %4.4f,',wcg);
    fprintf(' TotalCompTime = %4.4f\n',tTotal);
end
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
function [LMI,eLMI] = evalLMI(G,NE,v,p,g,X11,P,Pdot)
% This function evaluates the finite-horizon gain condition on a time grid
% for given values of the IQC, storage function, and gain. **The function
% assumes G is LTV and Delta is a SISO, unit norm bounded, LTI uncertainty.
%
% NOTE - All data (G,P,Pdot) must have the same time grid t

% Form Extended System
[AAe,BBe,CCe,DDe,~,Nz,~,Nd] = ltvutil.extsystem(G,v,p,NE);
n = size(AAe,1);
CC1e = CCe(1:Nz,:);
DD1e = DDe(1:Nz,:);
CC2e = CCe(Nz+1:end,:); % Should be zero for pure Euclidean Cost
DD2e = DDe(Nz+1:end,:); % Should be zero for pure Euclidean Cost

% Evaluate LMI
P.InterpolationMethod = 'Linear';
Pdot.InterpolationMethod = 'Linear';

X = blkdiag(X11,-X11);
gsq = g^2;
LMI = [AAe'*P+P*AAe+Pdot P*BBe; BBe'*P zeros(Nd+1)] ...
    + blkdiag( zeros(n+1), -gsq*eye(Nd) ) + [CC1e DD1e]'*X*[CC1e DD1e] ...
    + [CC2e DD2e]'*[CC2e DD2e];

% Evaluate max real eigenvalue
eLMI = max(real( eig(LMI) ));
end