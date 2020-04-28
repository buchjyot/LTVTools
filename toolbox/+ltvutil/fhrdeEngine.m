function [g,info] = fhrdeEngine(G,Delta,LMIinfo,Opt,NE)
% This function computes the finite-horizon gain for an uncertain system
% Fu(G,Delta) using the IQC LMI condition. **The function assumes G is LTI
% and Delta is a SISO, unit norm bounded, LTI uncertainty.
%
% Inputs
%   G - Nominal LTI system
%   (v,p) - The IQC multiplier is blkdiag(Psiv'*X11*Psiv,-Psiv'*X11*Psiv)
%           where Psiv:=[1; 1/(s-p); ...; 1/(s-p)^v] and X11>=0.
%   Tf - Finite horizon is [0,Tf]
%   X11 - IQC matrix
%   gBnds - 1-by-2 vector specifying initial bisection lower
%      bound gBnds(1) and upper bound gBnds(2).
%
% Outputs
%   g - Gain Upper bound
%   info - Data structure with information on RDE solutions

%% Input Processing
nin = nargin;
nout = nargout;
if nin<5
    NE = 0;
end
gLow = Opt.Bounds(1);
gUpp = Opt.Bounds(2);

%% Read Uncertainty Information
[Nv,Nw] = size(Delta);
[Type,Psi] = ltvutil.iqcEngine(Delta);
X11 = LMIinfo.X11;
if isequal(Type,2)
    % For time-varying X11 we need interpolation
    X11 = evalt(tvmat(X11),G.Time);
end
X = blkdiag(X11,-X11);

%% Options
Opt.OdeOptions = odeset('Events',@LOCALevents);
cdreOpt.OdeSolver = Opt.OdeSolver;
cdreOpt.OdeOptions = Opt.OdeOptions;
DispFlag = isequal(Opt.Display,'on');
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
RDEcnt = 0;
[T0,Tf] = getHorizon(G);

%% Form Extended System
[AAe,BBe,CCe,DDe,CEeT,Nz,~,Nd] = ltvutil.extsystem(G,Psi,Nv,Nw,NE);
n = size(AAe,1);
CC1e = CCe(1:Nz,1:n);
CC2e = CCe(Nz+1:end,1:n);
DD1e = DDe(1:Nz,:);
DD2e = DDe(Nz+1:end,:);

%% Build Cost Function Matrices
% Note:  R(t,g) = R0(t) - g^2 R1(t)
A = AAe;
B = BBe;
E = [];
Q = CC1e'*X*CC1e + CC2e'*CC2e;
S = CC1e'*X*DD1e + CC2e'*DD2e;
R0 = DD1e'*X*DD1e + DD2e'*DD2e;
R1 = blkdiag(0,eye(Nd));
F = CEeT'*CEeT;

%% Lower Bound Phase
%  Require R(t,g)=R0(t)-g^2*R1(t)<0 for all t in [0,T]
% XXX This assumes the time grid is sufficiently fine.
if gLow==0 % isempty(gLow)
    Nt = 1e3;
    t = linspace(T0,Tf,Nt);
    gLow = zeros(Nt,1);
    for i=1:Nt
        ti = t(i);
        R0i = tvsubs(R0,ti);
        R1i = R1;   % XXX: tvsubs(R1,ti);
        ev = eig(R0i,R1i);
        ev = ev( isfinite(ev) & ev>=0 );
        if isempty(ev)
            if DispFlag
                fprintf('\n R(t,g)<0 is not satisfied at t = %4.3f\n',ti);
            end
            return
        end
        gLow(i) = max(ev);
    end
    gLow = sqrt(max(gLow));
end
if DispFlag
    fprintf('\n Lower Bound = %4.3f',gLow);
end

%% Upper Bound Phase
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];
if isfinite(gUpp)
    % User specified a (finite) upper bound.
    % Verify (or disprove) this upper bound.
    gTry = gUpp;
    R = R0-gTry^2*R1;
    if nout==1
        P = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
        Pdot = []; sol = [];
    else
        [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
    end
    RDEcnt = RDEcnt + 1;
    
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
    % gFac = 10;
    gFac = 2;
    gUpp = gFac*gLow+1;  % XXX Better choice?
    
    cnt = 0;
    cntMax = 10;   % 8;
    haveUpper = false;
    while ~haveUpper && cnt<cntMax
        % Pick gamma
        cnt = cnt+1;
        gTry = gUpp;
        
        % Solve LTV Riccati Equation
        R = R0-gTry^2*R1;
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
        end
        RDEcnt = RDEcnt + 1;
        
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

%% Bisection Phase
if isfinite(gUpp)
    if DispFlag
        fprintf('\n');
    end
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        gTry = (gUpp+gLow)/2;
        
        % Solve LTV Riccati Equation
        R = R0-gTry^2*R1;
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],cdreOpt);
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
                fprintf(' gTry = %4.3f \t tConv= %4.3f  (Y) \t gLow = %4.3f \t gUpp = %4.3f\n', ...
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
end