function [K,CL,g,info] = tvhinfric(G,Ny,Nu,NE,gTry,Opt)
%% TVHINFSYN Engine
% RDE based Hinf Measurement Feedback Synthesis

% Number of outputs
nout = nargout;

% Read Options
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
DispFlag = isequal(Opt.Display,'on');
cdreOpt.OdeSolver = Opt.OdeSolver;
cdreOpt.OdeOptions = Opt.OdeOptions;
gLow = Opt.Bounds(1);
gUpp = Opt.Bounds(2);

% Horizon
[T0,Tf] = getHorizon(G);

%% Build Cost Function Matrices
% Store Original System
Gnom = G;

% Discriptor systems are not supported
E = []; %#ok<NASGU>
S = [];

% Extract syndata
[~,~,~,~,~,~,D12,D21,D22,CTf,DTf,~,GL2] = syndata(G,Ny,Nu,NE);
[Nx,~,Nd,NL2] = syniodim(G,Ny,Nu,NE);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE gain.');
end

% Assumptions/RankConditions checking
[r1,r2] = ltvutil.hinfrc(D12,D21,Ny,Nu);

% Normalize D12 and D21 matrices of GL2
[Gscl,r12,r21] = ltvutil.ioscale(GL2,Ny,Nu);

% Extract Scalled Plant Data
[A,B1,B2,C1,C2,D11,D12,D21] = syndata(Gscl,Ny,Nu);
B = [B1,B2];
C = [C1;C2];
D1o = [D11 D12];
Do1 = [D11;D21];

% Generic Cost Matrices
% CDRE of X
Rx = @(GAM) D1o'*D1o - diag([GAM^2*ones(1,Nd) zeros(1,Nu)]);
Qx = @(GAM) C1'*(eye(NL2)-D1o*(Rx(GAM)\D1o'))*C1;
Ax = @(GAM) A - B*(Rx(GAM)\D1o')*C1;
Fx = @(GAM) CTf'*(eye(NE)-DTf*(tvsubs(Rx(GAM),Tf)\DTf'))*CTf;

% CDRE of Y
Ry = @(GAM) Do1*Do1' - diag([GAM^2*ones(1,NL2) zeros(1,Ny)]);
Ay = @(GAM) A - (B1*Do1'*(Ry(GAM)\C));
Qy = @(GAM) B1*(eye(Nd)-Do1'*(Ry(GAM)\Do1))*B1';
Fy = zeros(Nx);

% Fy = zeros(Nx) implies M = inv(Fy) is Inf Fy = 0, means you know the
% initial states certainly, whereas if Fy > 0 then it means there could be
% uncertainity in the initial states, Fy = Inf means your initial state is
% not known completely which is an unlikely case.

% Wrap in to function handle
COSTfh = @(GAM) deal(Qx(GAM),Rx(GAM),S,Fx(GAM),Qy(GAM),Ry(GAM),S,Fy);

% Set RDE Count to 0 (RDECnt=1 includes X and Y both Riccati Equations)
RDEcnt = 0;

% Begin Timing
t0 = tic;

%% Lower Bound Phase

% Check the Feedthrough D11 condition for existance of a controller
% Ref. Theorem 17.1 (a)(i) from Zhou, K., Doyle, J. C., & Glover, K. (1996).
% Robust and optimal control (Vol. 40, p. 146). New Jersey: Prentice hall.
if any(D11.Data(:))
    % Partition D matrix so that D1111 is Np by Nq
    Np = NL2-r1;
    Nq = Nd-r2;
    D1111 = D11(1:Np,1:Nq);
    D1112 = D11(1:Np,Nq+1:Nd);
    D1121 = D11(Np+1:NL2,1:Nq);
    sig1 = svd([D1111 D1112]);
    sig2 = svd([D1111' D1121']);
    gLow = max([tvmax(sig1),tvmax(sig2),gLow]);
end

% Return if gTry is provided such that its less than lower bound
if ~isempty(gTry)
    if gTry < gLow
        % i.e. such gTry is not achievable
        return;
    end
end

if DispFlag && isempty(gTry)
    fprintf(' Lower Bound = %4.3f\n',gLow);
end

%% Upper Bound Phase

% Consider H2 Synthesis
% XXX:

% NOTE: Following code is commented out because it makes the entire
% function slow for NE = 0. We should also check that feedthough matrix is
% 0 for this synthesis problem.

% if isequal(NE,0)
%     warning('off','ltvtools:tvh2norm:noNL2erofeedthrough');
%     [~,CLH2]= tvh2syn(Gbfsc,Ny,Nu,Opt);
%     warning('on','ltvtools:tvh2norm:noNL2erofeedthrough');
%     gH2Hinf = tvnorm(CLH2);
%     gUpp = min([gH2Hinf(2),gUpp]);
% end

if ~isfinite(gUpp)
    % Instead of Inf start with some high finite value
    gUpp = 1e6;
end

if DispFlag && isempty(gTry)
    fprintf(' Upper Bound = %4.3f\n',gUpp);
end

%% GTRY
% Solve CDRE if gTry is not empty i.e. user specified gTry
if ~isempty(gTry)
    [Xinf,Yinf,Xdot,Ydot,Xsol,Ysol] = localSolveCDRE(gTry,Ax,B,Ay,C,COSTfh,T0,Tf,nout,cdreOpt);
    RDEcnt = RDEcnt + 1;
    
    % Store Info
    info.Try.X = Xinf;info.Try.Xdot = Xdot;info.Try.Xsol = Xsol;
    info.Try.Y = Yinf;info.Try.Ydot = Ydot;info.Try.Ysol = Ysol;
    
    % Return [] if not controller is not admissible
    if ~localCheckConvergence(gTry,Xinf,Yinf,COSTfh,T0,Tf)
        return;
    else
        info.Try.Gain = gTry; g = gTry;
    end
end

%% Bisection Phase
% Gamma Iterations if gTry is empty
if isempty(gTry)
    
    % Display
    if DispFlag
        fprintf(' ### Starting Bisection Phase:\n');
    end
    
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        % Bisect Gain
        gTry = (gLow+gUpp)/2;
        
        % Solve CDREs
        [X,Y,Xdot,Ydot,Xsol,Ysol] = localSolveCDRE(gTry,Ax,B,Ay,C,COSTfh,T0,Tf,nout,cdreOpt);
        RDEcnt = RDEcnt + 1;
        
        % Check convergence of X and Y
        if ~localCheckConvergence(gTry,X,Y,COSTfh,T0,Tf)
            % Either X or Y did not converge or coupling condition did not satisfy
            gLow = gTry;
            
            % Store Info
            info.Lower.Gain = gLow;
            info.Lower.X = X; info.Lower.Xdot = Xdot; info.Lower.Xsol = Xsol;
            info.Lower.Y = Y; info.Lower.Ydot = Ydot; info.Lower.Ysol = Ysol;
            
            if DispFlag
                fprintf(' gTry = %4.3f \t tConvX = %4.3f  \t tConvY = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,X.Time(1),Y.Time(1),gLow,gUpp);
            end
        else
            % X and Y converged and coupling condition met
            gUpp = gTry;
            
            % Store Info
            info.Upper.Gain = gUpp;
            info.Upper.X = X; info.Upper.Xdot = Xdot; info.Upper.Xsol = Xsol;
            info.Upper.Y = Y; info.Upper.Ydot = Ydot; info.Upper.Ysol = Ysol;
            
            if DispFlag
                fprintf(' gTry = %4.3f \t tConvX = %4.3f  \t tConvY = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,X.Time(1),Y.Time(1),gLow,gUpp);
            end
        end
    end
    if isfield(info,'Upper')
        g = gUpp;Xinf = info.Upper.X;Yinf = info.Upper.Y;
    else
        error('Did not find an upper bound.')
    end
    
    if isinf(g)
        return
    elseif DispFlag
        fprintf(' Final Bounds: [%.4f,%.4f]\n',gLow,gUpp);
    end
end

%% Design Hinf (sub) optimal measurement feedback central controller

% XXX You may formulate a secondary optimization problem which searches for
% the better controller that meets the objective

switch Opt.Method
    case 'BackOff'
        % Change the method to MinGamma so it can return the controller
        % with gBackOff performance
        OptBackOff = Opt;OptBackOff.Method = 'MinGamma';
        gBackOffTry = Opt.BackOffFactor*g;
        [K,CL,g,BackOffInfo] = tvhinfsyn(Gnom,Ny,Nu,NE,gBackOffTry,OptBackOff);
        info.BackOff = BackOffInfo;
        
    case 'MinGamma'
        % Evaluate Xinf and Yinf on the same time as the plant G time (legacy
        % approach)
        % Xinf = evalt(Xinf,G.Time);
        % Yinf = evalt(Yinf,G.Time);
        
        % Use Union of time grids
        Tgrid = union(union(Gscl.Time,Xinf.Time),Yinf.Time);
        
        % Evaluate on same time grid
        [RX,RY] = LOCALevalg(Rx,Ry,g);
        [Xinf,Yinf,RX,RY,r12,r21,Gscl] = evalt(Xinf,Yinf,RX,RY,r12,r21,Gscl,Tgrid);
        
        % Extract Plant Data on Tgrid
        [A,B1,B2,C1,C2,D11,D12,D21] = syndata(Gscl,Ny,Nu);
        B = [B1,B2];
        C = [C1;C2];
        D1o = [D11 D12];
        Do1 = [D11;D21];
        
        % Controller and Estimator Gains
        F = -RX\(D1o'*C1 + B'*Xinf);
        L = -(B1*Do1' + Yinf*C')/RY;
        F1inf = F(1:Nd,:);
        F2inf = F(Nd+1:Nu+Nd,:);
        L1inf = L(:,1:NL2);
        L2inf = L(:,NL2+1:Ny+NL2);
        
        % Partition D matrix so that D1111 is Np by Nq
        Np = NL2-r1;
        Nq = Nd-r2;
        D1111 = D11(1:Np,1:Nq);
        D1112 = D11(1:Np,Nq+1:Nd);
        D1121 = D11(Np+1:NL2,1:Nq);
        D1122 = D11(Np+1:NL2,Nq+1:Nd);
        
        % Partition Controller and Estimator gains
        % F11inf = F1inf(1:Nq,:); L11inf = L1inf(:,1:Np);
        F12inf = F1inf(Nq+1:end,:);
        L12inf = L1inf(:,Np+1:end);
        
        % TVSS matrices
        Nd11211 = size(D1121,1);
        Nd11122 = size(D1112,2);
        D11hat = -D1121*D1111'*((g^2*eye(Np)-D1111*D1111')\D1112) - D1122;
        D12hat = chol(eye(Nd11211)-D1121*((g^2*eye(Nq)-D1111'*D1111)\D1121'),'lower');
        D21hat = chol(eye(Nd11122)-D1112'*((g^2*eye(Np)-D1111*D1111')\D1112));
        Zinf = eye(Nx)-g^-2*Yinf*Xinf;
        
        % Book Results
        C2hat = -D21hat*(C2 + F12inf);
        C1hat = F2inf + D11hat*(D21hat\C2hat);
        B2hat = Zinf\(B2 + L12inf)*D12hat;
        B1hat = -Zinf\L2inf + B2hat*(D12hat\D11hat);
        Ahat = A + B*F + B1hat*(D21hat\C2hat);
        
        % All Controller Solutions
        AK = Ahat;
        BK = [B1hat B2hat];
        CK = [C1hat;C2hat];
        DK = [D11hat D12hat;D21hat zeros(size(D21hat,1),size(D12hat,2))];
        Minf = tvss(AK,BK,CK,DK);
        info.AS = Minf;
        
        % Get K from all Controllers
        mU = size(Minf,1) - Ny;
        nU = size(Minf,2) - Nu;
        K = tvss(lft(Minf.Data,zeros(mU,nU)),Minf.Time,Minf.InterpolationMethod);
        
        % Back-substitute scaling factors r12 and r21
        K = r12*K*r21;
        
        % If D22 feedthrough term exists then based on Page 454: Zhou, K., Doyle,
        % J., & Glover, K. (1996). Robust and optimal control.
        if any(D22.Data(:))
            K = feedback(K,evalt(D22,K.Time));
        end
        
        % Wrap around the LFT interconnection
        % XXX lft for tvss is slower, use lft of data instead and then wrap in
        % to tvss object
        Gnom = evalt(Gnom,K.Time);
        CL = tvss(lft(Gnom.Data,K.Data),K.Time,K.InterpolationMethod);
end

% Stop Timing
info.TotalTime = toc(t0);
info.RDEcnt = RDEcnt;
end

function [X,Y,Xdot,Ydot,Xsol,Ysol,RX,RY] = localSolveCDRE(gTry,Ax,B,Ay,C,COSTfh,T0,Tf,nout,cdreOpt)
%% LOCAL function to solve CDREs localSolveCDRE(gTry,Ax,B,C,COSTfh,E,T0,Tf,nout,cdreOpt);
% CDRE in X solved backwards in time (related to FI problem)
% CDRE in Y solved forward in time (related to Hinf estimation problem)

% Initialize Output
Xdot = [];
Xsol = [];
Ydot = [];
Ysol = [];

% Solve 2 LTV Riccati Equation
[QX,RX,SX,FX,QY,RY,SY,FY] = COSTfh(gTry);
[At,Ab] = LOCALevalg(Ax,Ay,gTry);
E = [];

if nout~=4
    % CDRE for X solve backwards in time
    X = cdre(At,B,QX,RX,SX,E,FX,[Tf T0],cdreOpt);
    % CDRE for Y solve forward in time
    Y = cdre(-Ab',C',-QY,-RY,SY,E,FY,[T0 Tf],cdreOpt);
else
    % CDRE for X solve backwards in time
    [X,~,Xdot,Xsol] = cdre(At,B,QX,RX,SX,E,FX,[Tf T0],cdreOpt);
    % CDRE for Y solve forward in time
    [Y,~,Ydot,Ysol] = cdre(-Ab',C',-QY,-RY,SY,E,FY,[T0 Tf],cdreOpt);
end
end

function out = localCheckConvergence(gTry,X,Y,COSTfh,T0,Tf)
%% Local function to check covergence of X and Y
% Returns true if converged & coupling condition met

% Check Convergence
XConvergence = ~(X.Time(1) > T0);
YConvergence = ~(Y.Time(end) < Tf);
out = XConvergence && YConvergence;

if out
    % Checking that X(0) < gamma^2*inv(FY)
    [~,~,~,~,~,~,~,FY] = COSTfh(gTry);
    if tvmin(tvmat(FY)) > 0
        % i.e. FY is Positive Definite Matrix and hence invertible
        X0 = tvsubs(X,X.Time(1));
        BoundryCondition = (tvmax(eig(X0*FY)) - gTry^2 < 0);
    else
        % As simple as X(0) < Inf
        BoundryCondition = out;
    end
    
    % Check Coupling Condition
    % rho(X(t)*Y(t)) < g^2
    if BoundryCondition
        X = evalt(X,Y.Time);
        rhoXY = tvmax(abs(eig(X*Y)));
        CouplingCondition = (rhoXY < gTry^2);
    else
        CouplingCondition = false;
    end
else
    return;
end

% Process Final Decision
out = out & BoundryCondition & CouplingCondition;
end

function varargout = LOCALevalg(varargin)
%% Evaluate the cost if it depends on gamma
% e.g. [Q,R,S,E,...] = evalg(Q,R,S,E,...,g)
nin = nargin;
g = varargin{end};
varargout = cell(1,nin-1);
for i = 1:nin-1
    if isa(varargin{i},'function_handle')
        fh = varargin{i};
        varargout{i} = fh(g);
    else
        varargout{i} = varargin{i};
    end
end
end