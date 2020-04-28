function [g,d,info] = tvnorm(G,varargin)
%% TVNORM Computes the time varying system norm
%
% [g,d,info] = tvnorm(G);
% Computes default 'L2toL2' norm with default options.
%
% [g,d,info] = tvnorm(G,Opt);
% Computes default 'L2toL2' with user specified options.
%
% [g,d,info] = tvnorm(G,NE);
% Computes tvnorm based on following cases with default options.
%
% [g,d,info] = tvnorm(G,NE,Opt);
% Computes tvnorm based on following cases with user specified Options.
%
% NOTE:
% Let G be a TVSS system with NY by NU IO dimentions
% NE is an optional argument with default set to 0.
% NE = NY means 'L2toE' Norm
% 0 < NE < NY then we have 'L2toE' and 'L2toL2' mixed cost terms.
%
% Output Arguments are worst-case gain (g) induced by L2 disturbance (d)
% and relevant Riccati solution info

%% Input Processing

% Check # of inputs
narginchk(1,3);
nin = nargin;
nout = nargout;

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
NY = size(G,1);

% Obtain default Opt
switch nin
    case 1
        Opt = tvnormOptions;
        NE = 0;
    case 2
        if isa(varargin{1},'tvnormOptions')
            Opt = varargin{1};
            NE = 0;
        else
            Opt = tvnormOptions;
            NE = varargin{1};
        end
    case 3
        NE = varargin{1};
        Opt = varargin{2};
end

% DispFlag
DispFlag = isequal(Opt.Display,'on');

% Update cdreOpt based on tvnormoptions
odeOpt = tvodeOptions('OdeOptions',Opt.OdeOptions,'OdeSolver',Opt.OdeSolver);

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% Counter for RDE bisections
RDEcnt = 0;

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY.');
end

%% Build Cost Function Matrices
[A,B] = ssdata(G);
[~,Nd] = size(B);

% We do not support discriptor systems currently
E = [];

% XXX (JB) The following should work when 1:NY-NE is empty by some size
% >> G(1:NY-NE,:)

% Extract L2 part
[~,~,CL2,DL2] = ssdata(tvss(G.Data(1:NY-NE,:),G.Time,G.InterpolationMethod));

% Extract Euclidean part
[~,~,CE,DE] = ssdata(tvss(G.Data(NY-NE+1:NY,:),G.Time,G.InterpolationMethod));
[CTf,DTf] = tvsubs(CE,DE,Tf);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm');
end

% Generic Cost Matrices
% XXX This set-up requires R(t) to be inverted at each time step.
% However, no inversion is actually required since we can set
% R=I and B'=B/gamma.  This should speed up the integration.
R0 = DL2'*DL2;
R1 = eye(Nd);
Rc = @(GAM) R0-(GAM^2)*R1;
Qc = CL2'*CL2;
Sc = CL2'*DL2;
Fc = CTf'*CTf;
COSTfh = @(GAM) deal(Qc,Rc(GAM),Sc,Fc);

%% Lower Bound Phase
% Require R(t,g)=R0(t)-g^2*R1(t)<0 for all t in [0,T]
% XXX This assumes the time grid is sufficiently fine.
evmax = tvmax(eig(R0));
gLow = sqrt(evmax);
gLow = max(gLow,Opt.Bounds(1));

% Display
if DispFlag
    fprintf(' Lower Bound = %4.3f\n',gLow);
end

%% Upper Bound Phase
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];
if isfinite(Opt.Bounds(2))
    % User specified a (finite) upper bound.
    % Verify (or disprove) this upper bound.
    gTry = Opt.Bounds(2);
    
    [Q,R,S,F] = COSTfh(gTry);
    if nout==1
        P = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
        Pdot = []; sol = [];
    else
        [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
    end
    
    % Check convergence of P
    if P.Time(1)>T0
        % P did not converge
        gUpp = inf;
        gLow = gTry;
        PLow = P; PdotLow = Pdot; solLow = sol;
    else
        % P converged
        gUpp = gTry;
        PUpp = P; PdotUpp = Pdot; solUpp = sol;
    end
else
    % User did not specify a (finite) upper bound.
    % Attempt to determine a finite upper bound on performance.
    gFac = 10;
    gUpp = gFac*gLow+1;  % XXX Better choice?
    
    cnt = 0;
    cntMax = 8;
    haveUpper = false;
    while ~haveUpper && cnt<cntMax
        % Pick gamma
        cnt = cnt+1;
        gTry = gUpp;
        
        % Solve LTV Riccati Equation
        [Q,R,S,F] = COSTfh(gTry);
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
        end
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gUpp = gFac*gUpp;
            gLow = gTry;
            PLow = P; PdotLow = Pdot; solLow = sol;
        else
            % P converged
            haveUpper = true;
            gUpp = gTry;
            PUpp = P; PdotUpp = Pdot; solUpp = sol;
        end
    end
    
    if ~haveUpper
        % Could not find a finite upper bound
        if DispFlag
            fprintf([' Could not find a finite upper bound.' ...
                ' Infeasible at gTry = %4.1f\n'],gUpp);
        end
        gUpp = inf;
    else
        if DispFlag
            fprintf(' Lower Bound = %4.3f \t Upper Bound = %4.3f\n',...
                gLow,gUpp);
        end
    end
end

if DispFlag
    fprintf(' Upper Bound = %4.3f\n',gUpp);
end

%% Bisection Phase
if DispFlag
    fprintf(' ### Starting Bisection Phase:\n');
end
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
if isfinite(gUpp)
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        gTry = (gUpp+gLow)/2;
        
        % Solve LTV Riccati Equation
        [Q,R,S,F] = COSTfh(gTry);
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],odeOpt);
        end
        RDEcnt = RDEcnt + 1;
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gLow = gTry;
            PLow = P; PdotLow = Pdot; solLow = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        else
            % P converged
            gUpp = gTry;
            PUpp = P; PdotUpp = Pdot; solUpp = sol;
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        end
        
    end
end

%% Store Final Result
g = [gLow,gUpp];

info.Lower.Gain = gLow;
info.Lower.P = PLow;
info.Lower.Pdot = PdotLow;
info.Lower.sol = solLow;

info.Upper.Gain = gUpp;
info.Upper.P = PUpp;
info.Upper.Pdot = PdotUpp;
info.Upper.sol = solUpp;

info.RDEcnt = RDEcnt;

%% Construct Worst-Case Input
% The construction is based on the two-point boundary value problem
% (TPBVP) related to the (Q,S,R,F) cost.  This particular implementation
% using the transformed Hamiltonian dynamics was developed by A. Iannelli
% and P. Seiler in September 2017.
d = [];
if nout>=2 && ~isempty(PLow)
    % Compute largest eigenvalue of PLow
    T = PLow.Time;
    [evecP,evalP] = eig( tvsubs(PLow,T(1)) );
    [emax,idx] = max( diag(evalP) );
    vmax = evecP(:,idx);
    
    %norm( tvsubs(PLow,T(1))*vmax - emax*vmax )  %norm(emax*vmax)]
    % Evaluate State/Cost Matrices on Same Time Grid as P
    [~,R,S,~] = COSTfh(gLow);
    [A,B,R,S] = evalt(A,B,tvmat(R),tvmat(S),T);
    
    % Simulate Transformed Hamiltonian system
    if isempty(S)
        M = R\(PLow*B)';
    else
        M = R\(PLow*B+S)';
    end
    H11 = A-B*M;
    odefh = @(t,E) tvsubs(H11,t)*E;
    E0 = vmax/emax;
    [tE,E] = ode45( odefh,T,E0);
    E = tvmat(E,tE);
    
    % Construct disturbance
    % Note: This calculation uses a non-convergent Riccati solution and
    % hence the disturbance starts at PLow.Time(1) > T0.
    d = -M*E;
    
    % Append zero input to disturbance on the window [PLow.Time(1) T0)
    dTime = d.Time;
    dData = d.Data;
    
    dTime = [T0; T0+0.999*(dTime(1)-T0); dTime];
    dData = cat(3,zeros(Nd,1,2), dData);
    d = tvmat(dData,dTime);
    
    % Normalize worst-case disturbance to have norm 1
    d = d/tvnorm(d);
end
end