function [K,CL,g,info] = tvhinfKFI(G,Nu,NE,gTry,Opt)
%% TVHINFFI ENGINE
% RDE based Hinf Full-Information Synthesis

% Number of outputs
nout = nargout;

% Init
K = [];
CL = [];
g = [];
info = [];

% Read Opts
gLow = Opt.Bounds(1);
gUpp = Opt.Bounds(2);
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
DispFlag = isequal(Opt.Display,'on');
cdreOpt.OdeSolver = Opt.OdeSolver;
cdreOpt.OdeOptions = Opt.OdeOptions;

% Time Horizon
[T0,Tf] = getHorizon(G);

%% Build Cost Function Matrices
% Store Original System
Gnom = G;

% Discriptor systems are not supported
E = [];

% Extract syndata
[~,~,~,~,~,~,D12,~,~,CTf,DTf,~,GL2] = syndata(G,0,Nu,NE);
[Nx,~,Nd,NL2] = syniodim(G,0,Nu,NE);

% Verify there is no feedthrough from d->e for L2toE
if any(DTf(:))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm');
end

% Assumptions/RankConditions checking
ltvutil.hinfrc(D12,[],0,Nu);

% Normalize D12 matrix of GL2
[Gscl,r12] = ltvutil.ioscale(GL2,0,Nu);

% Extract Scalled Plant Data
[A,B1,B2,C1,~,D11,D12] = syndata(Gscl,0,Nu);
B = [B1,B2];
D1o = [D11,D12];

% Cost Matrices
Sc = [];
Rc = @(GAM) D1o'*D1o - diag([GAM^2*ones(1,Nd) zeros(1,Nu)]);
Qc = @(GAM) C1'*(eye(NL2)-D1o*(Rc(GAM)\D1o'))*C1;
Ac = @(GAM) A - B*(Rc(GAM)\D1o')*C1;
Fc = @(GAM) CTf'*(eye(NE)-DTf*(tvsubs(Rc(GAM),Tf)\DTf'))*CTf;
COSTfh = @(GAM) deal(Qc(GAM),Rc(GAM),Sc,Fc(GAM));

% Set RDE Count to 0
RDEcnt = 0;

% Start Timing
t0 = tic;

%% Lower Bound Phase
% XXX Systematic way to verify (or disprove) this lower bound.

% Check the Feedthrough D11 condition for existance of a controller
% Ref. Theorem 17.6 (a) from Zhou, K., Doyle, J. C., & Glover, K. (1996).
% Robust and optimal control (Vol. 40, p. 146). New Jersey: Prentice hall.
if any(D11.Data(:))
    I = eye(size(D12,2));
    Dp = [I;zeros(size(I))];
    sig = svd(Dp'*D11);
    gLow = max([tvmax(sig),gLow]);
end

% Return if gTry is provided such that its less than lower bound
if ~isempty(gTry)
    if gTry < gLow
        % i.e. such gTry is not achievable
        return;
    end
end

if DispFlag && isempty(gTry)
    fprintf(' Lower Bound = %4.3f',gLow);
end

%% Upper Bound Phase
if isempty(gTry)
    if isfinite(gUpp)
        % User specified a (finite) upper bound.
        % Verify (or disprove) this upper bound.
        gTryUBP = gUpp;
        
        % Solve LTV Riccati Equation
        P = LOCALSolveCDRE(gTryUBP,Ac,B,COSTfh,E,T0,Tf,nout,cdreOpt);
        
        % Check convergence of P
        if P.Time(1)>T0
            % P did not converge
            gUpp = inf; gLow = gTryUBP;
        else
            % P converged
            gUpp = gTryUBP;
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
            gTryUBP = gUpp;
            
            % Solve LTV Riccati Equation
            P = LOCALSolveCDRE(gTryUBP,Ac,B,COSTfh,E,T0,Tf,nout,cdreOpt);
            
            % Check convergence of P
            if P.Time(1)>T0
                % P did not converge
                gUpp = gFac*gUpp;
                gLow = gTryUBP;
            else
                % P converged
                haveUpper = true;
                gUpp = gTryUBP;
            end
        end
        
        if ~haveUpper
            % Could not find a finite upper bound
            if DispFlag
                fprintf(['\n Could not find a finite upper bound.' ...
                    ' Infeasible at gTry = %4.1f\n'],gUpp);
            end
            gUpp = inf;
        else
            if DispFlag
                fprintf('\n Lower Bound = %4.3f \t Upper Bound = %4.3f',...
                    gLow,gUpp);
            end
        end
    end
    
    if DispFlag
        fprintf('\n Upper Bound = %4.3f \n',gUpp);
    end
end

%% GTRY
% Solve CDRE if gTry is not empty i.e. user specified gTry
if ~isempty(gTry)
    [P,Pdot,sol] = LOCALSolveCDRE(gTry,Ac,B,COSTfh,E,T0,Tf,nout,cdreOpt);
    info.Try.P = P;info.Try.Pdot = Pdot;info.Try.sol = sol;
    
    % Return [] if not controller is not admissible
    if P.Time(1) > T0
        return;
    else
        info.Try.Gain = gTry; g = gTry;
    end
end

%% Bisection Phase
% Gamma Iterations if gTry is empty
if isfinite(gUpp) && isempty(gTry)
    
    % Display
    if DispFlag
        fprintf(' ### Starting Bisection Phase:\n');
    end
    
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        % Bisect Gain
        gTry = (gLow+gUpp)/2;
        
        % Solve LTV Riccati Equation
        [P,Pdot,sol] = LOCALSolveCDRE(gTry,Ac,B,COSTfh,E,T0,Tf,nout,cdreOpt);
        RDEcnt = RDEcnt + 1;
        
        % Update gLow and gUpp
        % Check convergence of P
        if P.Time(1) > T0
            % P did not converge
            gLow = gTry;
            
            % Update Info
            info.Lower.P = P; info.Lower.Pdot = Pdot; info.Lower.sol = sol;
            info.Lower.Gain = gLow;
            
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (N) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        else
            % P converged
            gUpp = gTry;
            
            % Update Info
            info.Upper.Gain = gUpp;info.Upper.P = P;info.Upper.Pdot = Pdot;
            info.Upper.sol = sol;
            
            if DispFlag
                fprintf(' gTry = %4.3f \t tConv = %4.3f (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
                    gTry,P.Time(1),gLow,gUpp);
            end
        end
    end
    if isfield(info,'Upper')
        g = info.Upper.Gain;P = info.Upper.P;
    else
        error('Did not find an upper bound.')
    end
    
    if isinf(g)
        return
    elseif DispFlag
        fprintf(' Final Bounds: [%.4f,%.4f]\n',gLow,gUpp);
    end
end

%% Design Hinf (sub) optimal full state feedback central controller

% XXX You may formulate a secondary optimization problem which searches for
% the better controller that meets the objective.
% Evaluate B2 on P time grid as Riccati Solution may be on finer time grid
% compared to system time grid (Legacy Code)
% [D1o,C1,A,B,R] = LOCALEvalt(D1o,C1,A,B,Rx(g),RP.Time);

switch Opt.Method
    case 'BackOff'
        % Change the method to MinGamma so it can return the controller
        % with gBackOff performance
        OptBackOff = Opt;OptBackOff.Method = 'MinGamma';
        gBackOffTry = Opt.BackOffFactor*g;
        [K,CL,g,BackOffInfo] = tvhinfsyn(Gnom,Ny,Nu,NE,gBackOffTry,OptBackOff);
        info.BackOff = BackOffInfo;
        
    case 'MinGamma'
        % Use Union of time grids
        Tgrid = union(Gscl.Time,P.Time);
        
        % Make all matrices on common time grid
        [~,R,~,~] = COSTfh(g);
        [P,Gscl,Gnom,R,r12] = evalt(P,Gscl,Gnom,tvmat(R),r12,Tgrid);
        
        % Extract Plant Data on Tgrid for controller reconstruction
        [~,B1,B2,C1,~,D11,D12] = syndata(Gscl,0,Nu);
        B = [B1 B2];
        D1o = [D11 D12];
        
        % Designed Hinf central state feedback controller Theorem 17.6 from
        % ZDG Book (NOTE: Put back the scaling factor r12)
        Kinf = -(R\(D1o'*C1 + B'*P));
        FT = D12'*D11;
        Kd = r12*(-FT);
        Kx = r12*[FT  eye(size(FT,1))]*Kinf;
        K = [Kx Kd];
        
        % Generate Extended System based off Gnom
        [Anom,Bnom,C1ext,D1ext] = ssdata(Gnom);
        Cext = [C1ext;eye(Nx);zeros(Nd,Nx)];
        Dext = [D1ext;zeros(Nx,Nd+Nu);[eye(Nd) zeros(Nd,Nu)]];
        Greg = tvss(Anom,Bnom,Cext,Dext);
        Greg = evalt(Greg,K.Time);
        
        % Form closed loop using LFT Interconnection
        CL = lft(Greg,K);
        info.Greg = Greg;
        info.Kx = Kx;
        info.Kd = Kd;
end

info.TotalTime = toc(t0);
info.RDEcnt = RDEcnt;
end

function [P,Pdot,Psol,R] = LOCALSolveCDRE(gTry,Ax,Bx,COSTfh,E,T0,Tf,nout,cdreOpt)
%% Local function to solve the CDRE
% Solving the following LTV Riccati equation backwards in time
% - Pdot = At'*P + P*At - P*(B2*B2' - g^-2*B1*B1')*P + Ct'*Ct, P(T) = F

% Initialize Outputs
Pdot = [];
Psol = [];

% Gamma specific cost matrices
[Q,R,S,F] = COSTfh(gTry);
[At,B] = LOCALevalg(Ax,Bx,gTry);

if nout~=4
    P = cdre(At,B,Q,R,S,E,F,[Tf T0],cdreOpt);
else
    [P,~,Pdot,Psol] = cdre(At,B,Q,R,S,E,F,[Tf T0],cdreOpt);
end
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