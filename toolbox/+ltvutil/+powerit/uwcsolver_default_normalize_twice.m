function [glb,wcInput,info] = uwcsolver_default_normalize_twice(mySystem,tUser,U1,x0,pSpec,Opt,systemType)
%% UWCSOLVER Routine
%
% This function assumes that uncertainty is wrapped on the top of the plant
% i.e. lft(Delta,mySystem);

% Requested number of output arguments
nout = nargout;

% Populate signal sizes
fnames = fieldnames(pSpec);
Nfield = length(fnames);
for k = 1:Nfield
    eval([fnames{k} '= pSpec.(fnames{k});']);
end

% We allow only square full block uncertainty
if ~isequal(Nv,Nw)
    error('Dynamic uncertainty must be square i.e. Nw = Nv.');
end

% No parametric uncertainty supported
if ~isequal(Np,0)
    error('Parametric uncertainties are not supported yet.');
end

% Only TVSS systems are supported
if ~isequal(systemType,'tvss')
    error('Only time-varying state-space systems are supported.')
end

% Identify Options
MaxIter          = Opt.MaxIter;
StopTol          = Opt.StopTol;
OdeSolver        = Opt.OdeSolver;
OdeOptions       = Opt.OdeOptions;
LinOpt           = Opt.LinOpt;
SimOpt           = Opt.SimOpt; %#ok<*NASGU>
Display          = Opt.Display;
StoreAllIter     = (Opt.StoreAllIter && nout>2);
NumPertRel       = LinOpt.NumericalPertRel;
tvopt            = tvodeOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions);

% Create a watch window
if isequal(Display,'plot')
    myWatcher = poweritWatcher;
else
    myWatcher = [];
end

%% Normalize Initial Input

% Uncertain input
w = U1(1:Nw,:);
w = w*(UncNormBnd/tvnorm(w));

% Uncertain output
u = U1(Nw+1:Nw+Nu,:);
u = u*(InputL2Norm/tvnorm(u));

% Precompute scallings for first iteration
normw = tvnorm(w);
normu = tvnorm(u);

% Combine input
U1 = [w;u];

%% PlaceHolders
% Plant Signals
thisU    = U1;
thisY    = [];
thisX    = [];

% Adjoint Signals
thisUa   = []; %#ok<*NASGU>
thisYa   = U1;
thisXa   = [];

% Change in input signal
deltaU   = [];

% Keep track of previous iteration signals
prevU    = [];
prevY    = [];
prevX    = [];

% Keep track of worst-case signals
wcU      = [];
wcY      = [];
wcX      = [];

% Keep track of current, previous and worst-case performance
thisPerf = [];
prevPerf = [];
wcPerf   = [];

% Preallocate memory
if StoreAllIter
    AllIter(MaxIter,1) = struct('FwdPerf',[],'U',[],'Y',[],'X',[],'AdjPerf',[],'Ua',[],'Ya',[],'Xa',[]);
end
DisplayString = @ltvutil.powerit.dispstr;

%% Setup

% Turn off all warnings
orig_warning_state = warning('off','all');

%% Power Iterations Main loop

% Display starting
[T0,Tf0] = getHorizon(U1);
if isequal(Display,'on')
    fprintf(' Starting power iterations...\n');
    fprintf(' Analysis horizon is fixed to [%.3f, %.3f] seconds.\n',T0,Tf0);
end

% Starting timing
t1 = tic;

% Start iterations
for i = 1:MaxIter
    
    %% Simulate the system
    % Solve 2*Nx linear ode equations
    thisU1 = [w;zeros(Nu,1)]; % Uncertain Inputs
    thisU2 = [zeros(Nw,1);u]; % Performance Inputs
    [thisY1,thisX1] = tvlsim(mySystem,thisU1,tUser,x0,tvopt); % Uncertain Input Response
    [thisY2,thisX2] = tvlsim(mySystem,thisU2,tUser,x0,tvopt); % Performance Input Response
    [wcg1,v,yL2,yE,thisY,thisX] = LOCALs4vecp(thisY1,thisY2,thisX1,thisX2,Nv,NL2,NE,Tf0);
    
    %% Linearize System
    switch systemType
        case 'tvss'
            % TVSS objects are already linear system, thus you need to do
            % this only the first time you run this loop.
            if isequal(i,1)
                [dfdx,dfdu,dgdx,dgdu] = ssdata(mySystem);
            end
    end
    
    %% Costate Dynamics and Boundary Condition
    DTf     = tvsubs(dgdu(Nv+NL2+1:end,:),Tf0);
    if any(DTf(:))
        error('Feedthrough term d->e must be zero for well-posed L2toE performance.');
    end
    Ga      = tvss(-dfdx',-dgdx(1:Nv+NL2,:)',dfdu',dgdu(1:Nv+NL2,:)');
    lamTf   = tvsubs(dgdx(Nv+NL2+1:end,:),Tf0)'*yE;
    
    %% Allignment Condition 1
    % Compute norms
    normv   = tvnorm(v);
    normyL2 = tvnorm(yL2);
    normyE  = norm(yE);
    normy   = sqrt(normyE^2 + normyL2^2);
    
    % Adjoint inputs
    thisUa1 = [v*(normw/normv);zeros(NL2,1)];
    thisUa2 = [zeros(Nv,1);yL2*(normu/normy)];
    thisUa  = thisUa1 + thisUa2;
    
    %% Simulate Costate Dynamics backwards in Time
    % Solve 2*Nx linear ode equations for adjoint
    [thisYa1,thisXa1] = tvlsim(Ga,thisUa1,flip(tUser),lamTf,tvopt);
    [thisYa2,thisXa2] = tvlsim(Ga,thisUa2,flip(tUser),lamTf,tvopt);
    [wcg2,w,u,~,thisYa,thisXa] = LOCALs4vecp(thisYa1,thisYa2,thisXa1,thisXa2,Nw,Nu,0,T0);
    
    %% Store This Iterations
    % Pete says that in most cases wcg1 and wcg2 should converge to the
    % same value,in any case, if wcgain1 and wcgain2 don't converge then we
    % need to think more carefully about how to interpret the result.
    
    thisPerf = wcg1;% max(wcg1,wcg2);
    if StoreAllIter
        % System performance and signals
        AllIter(i).FwdPerf  = thisPerf;
        AllIter(i).U        = thisU;
        AllIter(i).X        = thisX;
        AllIter(i).Y        = thisY;
        
        % Adjoint performance and signals
        AllIter(i).FwdPerf  = wcg2;
        AllIter(i).Ua       = thisUa;
        AllIter(i).Xa       = thisXa;
        AllIter(i).Ya       = thisYa;
    end
    
    %% Stopping Condition
    if (i > 1) && (i < MaxIter)
        % Keep track of maximum gain and related signals
        [gMax,id] = max([prevPerf,thisPerf,wcPerf]);
        switch id
            case 1
                % Previous iteration was maximum
                wcU     = prevU;
                wcX     = prevX;
                wcY     = prevY;
                wcPerf  = prevPerf;
                
            case 2
                % This iteration is maximum
                wcU     = thisU;
                wcX     = thisX;
                wcY     = thisY;
                wcPerf  = thisPerf;
                
            case 3
                % No change to worst-case iteration
        end
        
        % Check if performance that we care is stationary
        if abs(thisPerf-prevPerf) <= StopTol*gMax
            % Check if input signal is stationary
            [StopTolSatisfied,deltaU] = ltvutil.powerit.stopcond(StopTol,thisU,prevU);
            if StopTolSatisfied
                break;
            end
        end
    end
    
    %% Display
    switch Display
        case 'on'
            fprintf(DisplayString(i,thisPerf,deltaU));
            fprintf(newline);
            
        case 'plot'
            myWatcher = updatewatch(myWatcher,thisU,DisplayString(i,thisPerf,deltaU));
            if hasterminated(myWatcher)
                break;
            end
    end
    
    %% Update This Iteration as Previous Iteration
    prevPerf = thisPerf;
    prevU    = thisU;
    prevX    = thisX;
    prevY    = thisY;
    
    %% Allignment Condition 2
    % Compute Norms
    normw = tvnorm(w);
    normu = tvnorm(u);
    
    % Input for the next iteration
    w = w*(normv/normw);
    u = u*(normy/normu);
    thisU = [w; u];
end

% Stop timing
tTotal = toc(t1);

%% Cleanup

% Turn on all warning
warning(orig_warning_state);

%% Terminate display
switch Display
    case 'on'
        % Terminate command window display
        if isequal(i,MaxIter) && ~StopTolSatisfied
            fprintf(' Maximum number of iterations reached.\n');
        elseif StopTolSatisfied
            fprintf(DisplayString(i,thisPerf,deltaU));fprintf(newline);
            fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
        end
        fprintf(' Final Results: NominalPerfLB = %4.4f, TotalCompTime = %4.4f\n',wcPerf,tTotal);
        
    case 'plot'
        % Shutdown the watch window's interactive features
        myWatcher = updatewatch(myWatcher,thisU,DisplayString(i,thisPerf,deltaU));
        closewatch(myWatcher);
end

%% Final Output
info = [];
if isequal(i,MaxIter) && ~StopTolSatisfied
    % In this case, algorithm did not converge, we should return the
    % worst-case iteration info.
    glb       = wcPerf;
    wcInput   = wcU;
    info.Xwc  = wcX;
    info.Ywc  = wcY;
    
elseif StopTolSatisfied
    % If stopping tolerance was satisfied then return the performance
    % achieved by stationary input
    glb       = thisPerf;
    wcInput   = thisU;
    info.Xwc  = thisX;
    info.Ywc  = thisY;
end

% Return general information
info.TotalTime   = tTotal;
info.TotalIter   = i;
if StoreAllIter
    info.AllIter = AllIter(1:i);
    if StopTolSatisfied
        info.WCIter = i;
    else
        FwdPerf      = [AllIter.FwdPerf];
        [~,Nwc]      = max(FwdPerf);
        info.WCIter  = Nwc;
    end
end
end

function [wcg,v,yL2,yE,thisY,thisX] = LOCALs4vecp(thisY1,thisY2,thisX1,thisX2,Nv,NL2,NE,Tterm)
%% LOCALs4vecp
% Subroutine for SKEWED power method.
% Find alpha (or 1/wcg) such that tvnorm of output is 1

% Evaluate on same time grid
tUnion = union(thisY1.Time,thisY2.Time);
[thisY1,thisY2,thisX1,thisX2] = evalt(thisY1,thisY2,thisX1,thisX2,tUnion);

% Construct a1,a2,b1,b2,c1,c2 such that
% v   = a1 + a2*sqrt(alpha)
% yL2 = sqrt(alpha)*(b1 + b2*sqrt(alpha))
% yE  = sqrt(alpha)*(c1(T) + c2(T)*sqrt(alpha));

a1 = thisY1(1:Nv,:);
b1 = thisY1(Nv+1:Nv+NL2,:);
c1 = tvsubs(thisY1(Nv+NL2+1:Nv+NL2+NE,:),Tterm);

a2 = thisY2(1:Nv,:);
b2 = thisY2(Nv+1:Nv+NL2,:);
c2 = tvsubs(thisY2(Nv+NL2+1:Nv+NL2+NE,:),Tterm);

% Compute uncertain channel tvnorms
na1 = tvnorm(a1);
na2 = tvnorm(a2);

% Compute L2 performance channel tvnorms
nb1 = tvnorm(b1);
nb2 = tvnorm(b2);

% Compute Euclidean performance channel norms
nc1 = norm(c1);
nc2 = norm(c2);

% Euclidean cross terms
c1c2 = c1'*c2;

% Compute L2 performance channel cross terms
b1b2 = trapz(reshapedata(b1'*b2),tUnion);

% Compute Uncertain channel cross terms
a1a2 = trapz(reshapedata(a1'*a2),tUnion);

% NOTE: First we solve the associated quartic polynomial for norm of output
% being 1. Due to numerical integration errors and time gridding this may
% not give exatly norm 1. If there are discripancies then we do bisection
% to resolve them with some numerical tolerance.

% Quartic poly in sqrt(alpha) for terminalPenalty + integralExpression - 1 = 0
betavec = [nc2^2+nb2^2, 2*(b1b2+c1c2), na2^2+nb1^2+nc1^2, 2*a1a2, na1^2-1];
tmp = ltvutil.roots(betavec);
ind = find( imag(tmp)==0 & tmp>0, 1 );
if ~isempty(ind)
    % Results based on quartic poly
    beta = tmp(ind(1));
    alpha = beta^2;
    
    % Worst-case gain
    wcg  = 1/alpha;
    
    % Construct outputs
    sa  = sqrt(alpha);
    v   = a1+a2*sa;
    yL2 = sa*(b1+b2*sa);
    yE  = sa*(c1+c2*sa);
    
    % Check if the output performance is unity
    tvn = sqrt(norm(yE)^2 + tvnorm([v;yL2])^2); % must be 1, can have numerical errors
    tvnTol = 1e-1;
    if abs(tvn-1) > tvnTol
        gLowIn = 0;
        if tvn > 1
            gUppIn = alpha;
        else
            gUppIn = alpha*1e3;
        end
        [sa,wcg,v,yL2,yE] = LOCALBisect(a1,a2,b1,b2,c1,c2,gLowIn,gUppIn);
    end
else
    [sa,wcg,v,yL2,yE] = LOCALBisect(a1,a2,b1,b2,c1,c2,0,1e6);
end

% Due to linearity perform super-position
thisX  = thisX1 + thisX2*sa;
OutScl = blkdiag(eye(Nv),sa*eye(NL2+NE));
thisY  = OutScl*(thisY1 + thisY2*sa);
end

function [sa,wcg,v,yL2,yE] = LOCALBisect(a1,a2,b1,b2,c1,c2,gLowIn,gUppIn)
%% LOCALBisect
gLow = max(0,gLowIn);
gUpp = min(gUppIn,1e6); % XXX Better choice?
BisectionTol = 1e-5;
cnt = 0;
while abs(gUpp-gLow) > BisectionTol
    % Pick gTry
    gTry = 0.5*(gUpp+gLow);
    sg = sqrt(gTry);
    
    % Construct outputs
    yE  = sg*(c1+c2*sg);
    yL2 = sg*(b1+b2*sg);
    v   = a1+a2*sg;
    
    % Bisect
    if sqrt(norm(yE)^2 + tvnorm([v;yL2])^2)<1
        gLow = gTry;
    else
        gUpp = gTry;
    end
    cnt = cnt + 1;
end
alpha = gTry;

% Worst-case gain
wcg  = 1/alpha;

% Construct outputs
sa  = sqrt(alpha);
v   = a1+a2*sa;
yL2 = sa*(b1+b2*sa);
yE  = sa*(c1+c2*sa);
end