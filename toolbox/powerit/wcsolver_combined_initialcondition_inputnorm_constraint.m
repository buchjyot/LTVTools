function [glb,wcInput,info] = wcsolver_combined_initialcondition_inputnorm_constraint(mySystem,tUser,U1,x0,pSpec,Opt,systemType)
%% WCSOLVER Routine

% Requested number of output arguments
nout = nargout;

% Populate signal sizes
fnames = fieldnames(pSpec);
Nfield = length(fnames);
for k = 1:Nfield
    eval([fnames{k} '= pSpec.(fnames{k});']);
end

% Identify Options
MaxIter          = Opt.MaxIter;
StopTol          = Opt.StopTol;
OdeSolver        = Opt.OdeSolver;
OdeOptions       = Opt.OdeOptions;
LinOpt           = Opt.LinOpt;
SimOpt           = Opt.SimOpt;
Display          = Opt.Display;
StoreAllIter     = (Opt.StoreAllIter && nout>2);
NumPertRel       = LinOpt.NumericalPertRel;
tvopt            = tvodeOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions);
StopTolSatisfied = false;

% Create a watch window
if isequal(Display,'plot')
    myWatcher = poweritWatcher;
else
    myWatcher = [];
end

% Initial Conditions Cost Matrix
E0 = InitialCondCostMat;
iE0 = inv(E0);

% Is initial condition fixed?
fixedIC = isequal(InitialConditions,'fixed');
if fixedIC
    DispString = @(in1,in2,in3,in4) dispstr(in1,in2,in3,[]);
else
    DispString = @(in1,in2,in3,in4) dispstr(in1,in2,in3,in4);
end

%% PlaceHolders
% Plant Signals
thisU    = U1;
thisY    = [];
thisX    = [];
thisX0   = x0;

% Adjoint Signals
thisUa   = []; %#ok<*NASGU>
thisYa   = [];
thisXa   = [];

% Change in input signal
deltaU   = [];
dX0      = [];

% Keep track of previous iteration signals
prevU    = [];
prevY    = [];
prevX    = [];
prevX0   = [];

% Keep track of worst-case signals
wcU      = [];
wcY      = [];
wcX      = [];
wcX0     = [];

% Keep track of current, previous and worst-case performance
thisPerf = [];
prevPerf = [];
wcPerf   = [];

% Preallocate memory
if StoreAllIter
    AllIter(MaxIter,1) = struct('FwdPerf',[],'U',[],'Y',[],'X',[],'AdjPerf',[],'Ua',[],'Ya',[],'Xa',[]);
end

%% Setup
switch systemType
    case 'model'
        % Simulink Model Setup: Configure model to have initial input
        % supplied externally through root level inports. Power iteration
        % code in the next subsection constructs the simInput variable in
        % the "base" workspace that is set as an external input to the
        % model inports.
        tempIn = createInputDataset(mySystem);
        set_param(mySystem,'LoadExternalInput','on');
        set_param(mySystem,'ExternalInput','simInput');
        
    case 'polysys'
        [~,~,AFunc,BFunc,CFunc,DFunc] = function_handle(mySystem);
        PartialsFunc = {AFunc,BFunc,CFunc,DFunc};
        
    case 'sfunc'
        PartialsFunc =  Opt.PartialsFunc;
end

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

% Start Iterations
for i = 1:MaxIter
    
    %% Simulate the system
    % Simulate state equation forward on a "Fixed" horizon
    switch systemType
        case 'tvss'
            [thisY,thisX] = tvlsim(mySystem,thisU,tUser,thisX0,tvopt);
            Tgrid1 = thisY.Time;
            
        case 'model'
            [tempInData,tempInTime] = reshapedata(thisU);
            tempIn{1} = timeseries(tempInData,tempInTime);
            assignin('base','simInput',tempIn);
            [Tgrid1,x,y] = sim(mySystem,tUser,SimOpt);
            [thisY,thisX] = double2tvmat(y,x,Tgrid1);
            
        case 'polysys'
            [tempInData,tempInTime] = reshapedata(thisU);
            [Tgrid1,x,y]  = sim(mySystem,tUser,thisX0,[tempInTime tempInData],OdeOptions,OdeSolver);
            [thisY,thisX] = double2tvmat(y,x,Tgrid1);
            
        case 'sfunc'
            [Tgrid1,thisX,thisY] = simsfunc(mySystem,tUser,thisU,thisX0,tvopt);
    end
    
    %% Compute Performance
    % Related signals
    yL2      = thisY(1:NL2,:);
    yE       = tvsubs(thisY(NL2+1:end,:),Tf0);
    
    % Compute norms
    normyL2  = tvnorm(yL2);
    normyE   = norm(yE);
    normy    = sqrt(normyE^2 + normyL2^2);
    nU       = tvnorm(thisU);
    nx0      = thisX0'*E0*thisX0;
    normu    = sqrt(nx0 + nU^2);
    
    % Compute induced gain
    FwdPerf  = normy/normu;
    
    %% Linearize System
    Nt = length(Tgrid1);
    switch systemType
        case 'tvss'
            % TVSS objects are already linear system, thus you need to do
            % this only the first time you run this loop.
            if isequal(i,1)
                [dfdx,dfdu,dgdx,dgdu] = ssdata(mySystem);
            end
            
        case 'model'
            % Compute costate dynamics using snapshot linearization about
            % the nominal trajectory thisU and thisX
            [linsys,op,lininfo]   = linearize(mySystem,Tgrid1,LinOpt); %#ok<ASGLU>
            Tgrid2                = linsys.SamplingGrid.Time;
            Glin                  = tvss(linsys,Tgrid2);
            [dfdx,dfdu,dgdx,dgdu] = ssdata(Glin);
            
            % Throw error if the simulation and linearization horizon(s)
            % are different from each other. This can happen if the
            % simulink model is configured to stop at some terminal event,
            % rather than a fixed horizon.
            Tf1 = Tgrid2(end);
            if ~isequal(Tf0,Tf1)
                error(['The horizon must be fixed. The nominal horizon was %.3f ',...
                    'seconds\nwhereas, the linearization was obtained on',...
                    'horizon %.3f seconds'],Tf0,Tf1);
            end
            
        case 'polysys'
            % Extract Data
            XiData = reshapedata(thisX);
            UiData = reshapedata(evalt(thisU,Tgrid1)); % Interpolated inputs
            
            % Memory allocation
            dfdx = zeros(Nx,Nx,Nt);
            dfdu = zeros(Nx,NIn,Nt);
            dgdx = zeros(NOut,Nx,Nt);
            dgdu = zeros(NOut,NIn,Nt);
            
            % Linearize
            for j = 1:Nt
                dfdx(:,:,j) = feval(PartialsFunc{1},Tgrid1(j),XiData(j,:)',UiData(j,:)');
                dfdu(:,:,j) = feval(PartialsFunc{2},Tgrid1(j),XiData(j,:)',UiData(j,:)');
                dgdx(:,:,j) = feval(PartialsFunc{3},Tgrid1(j),XiData(j,:)',UiData(j,:)');
                dgdu(:,:,j) = feval(PartialsFunc{4},Tgrid1(j),XiData(j,:)',UiData(j,:)');
            end
            [dfdx,dgdx,dfdu,dgdu] = double2tvmat(dfdx,dgdx,dfdu,dgdu,Tgrid1);
            
        case 'sfunc'
            % Extract Data
            XiData = reshapedata(thisX);
            UiData = reshapedata(evalt(thisU,Tgrid1)); % Interpolated inputs
            
            % Memory allocation
            dfdx = zeros(Nx,Nx,Nt);
            dfdu = zeros(Nx,NIn,Nt);
            dgdx = zeros(NOut,Nx,Nt);
            dgdu = zeros(NOut,NIn,Nt);
            
            % Linearize
            for j = 1:Nt
                linParameters = [ NumPertRel; Tgrid1(j); 0 ];
                if isempty(PartialsFunc)
                    [dfdx(:,:,j),dfdu(:,:,j),dgdx(:,:,j),dgdu(:,:,j)] = linsfunc(mySystem,XiData(j,:),UiData(j,:),linParameters);
                else
                    % Explicitly calculate partials
                    [dfdx(:,:,j),dfdu(:,:,j),dgdx(:,:,j),dgdu(:,:,j)] = feval(PartialsFunc,Tgrid1(j),XiData(j,:),UiData(j,:));
                end
            end
            [dfdx,dgdx,dfdu,dgdu] = double2tvmat(dfdx,dgdx,dfdu,dgdu,Tgrid1);
            
    end
    
    %% Costate Dynamics and Boundary Condition
    DTf     = tvsubs(dgdu(NL2+1:end,:),Tf0);
    if any(DTf(:))
        error('Feedthrough term d->e must be zero for well-posed L2toE performance.');
    end
    Ga      = tvss(-dfdx',-dgdx(1:NL2,:)',dfdu',dgdu(1:NL2,:)');
    lamTf   = tvsubs(dgdx(NL2+1:end,:),Tf0)'*yE;
    thisUa  = yL2;
    
    %% Simulate Costate Dynamics backwards in Time
    [thisYa,thisXa] = tvlsim(Ga,thisUa,flip(tUser),lamTf,tvopt);
    
    %% Compute Adjoint Performance
    % Compute input norm
    nlamTf = norm(lamTf);
    nUaL2  = tvnorm(thisUa);
    normUa = sqrt(nlamTf^2 + nUaL2^2);
    
    % Energy in initial condition
    lam0   = tvsubs(thisXa,T0);
    nlam0  = lam0'*iE0*lam0; %#ok<MINV>
    
    % Compute output signal norm
    nYa    = tvnorm(thisYa);
    normYa = sqrt(nlam0 + nYa^2);
    
    % Compute induced gain
    AdjPerf = normYa/normUa;
    
    %% Store This Iterations
    thisPerf = FwdPerf; % max(FwdPerf,AdjPerf);
    if StoreAllIter
        % System performance and signals
        AllIter(i).FwdPerf = FwdPerf;
        AllIter(i).U       = thisU;
        AllIter(i).X       = thisX;
        AllIter(i).Y       = thisY;
        
        % Adjoint performance and signals
        AllIter(i).AdjPerf = AdjPerf;
        AllIter(i).Ua      = thisUa;
        AllIter(i).Xa      = thisXa;
        AllIter(i).Ya      = thisYa;
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
                wcX0    = prevX0;
                wcPerf  = prevPerf;
                
            case 2
                % This iteration is maximum
                wcU     = thisU;
                wcX     = thisX;
                wcY     = thisY;
                wcX0    = thisX0;
                wcPerf  = thisPerf;
                
            case 3
                % No change to worst-case iteration
        end
        
        % Check if performance that we care is stationary
        if abs(thisPerf-prevPerf) <= StopTol*gMax
            % Check if input signal is stationary
            [StopTolSatisfied,deltaU,dX0] = stopcond(StopTol,thisU,prevU,thisX0,prevX0);
            if StopTolSatisfied
                break;
            end
        end
    end
    
    %% Display
    switch Display
        case 'on'
            fprintf(DispString(i,thisPerf,deltaU,dX0));
            fprintf(newline);
            
        case 'plot'
            myWatcher = updatewatch(myWatcher,thisU,DispString(i,thisPerf,deltaU,dX0));
            if hasterminated(myWatcher)
                break;
            end
    end
    
    %% Update This Iteration as Previous Iteration
    prevU    = thisU;
    prevX    = thisX;
    prevY    = thisY;
    prevX0   = thisX0;
    prevPerf = thisPerf;
    
    %% Allignment Condition
    % Compute next iteration system input
    thisU = thisYa*(InputL2Norm/normYa);
    
    % Update initial condition and normalize
    if ~fixedIC
        thisX0 = iE0*lam0/normYa; %#ok<MINV>
    end
end

% Stop timing
tTotal = toc(t1);

%% Cleanup
% Simulink Model Cleanup
if isequal(systemType,'model')
    % Turn off root level input
    set_param(mySystem,'LoadExternalInput','off');
    
    % Clear temporary base workspace variable
    evalin('base','clear simInput');
end

% Turn on all warning
warning(orig_warning_state);

%% Terminate display
switch Display
    case 'on'
        % Terminate command window display
        if isequal(i,MaxIter) && ~StopTolSatisfied
            fprintf(' Maximum number of iterations reached.\n');
        elseif StopTolSatisfied
            fprintf(DispString(i,thisPerf,deltaU,dX0));fprintf(newline);
            fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
        end
        fprintf(' Final Results: NominalPerfLB = %4.4f, TotalCompTime = %4.4f\n',wcPerf,tTotal);
        
    case 'plot'
        % Shutdown the watch window's interactive features
        myWatcher = updatewatch(myWatcher,thisU,DispString(i,thisPerf,deltaU,dX0));
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