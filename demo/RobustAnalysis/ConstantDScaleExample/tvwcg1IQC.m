function [g,d,info] = tvwcg1IQC(G,IQC,Tspan,varargin)
%% XXX Add Documentation
%   [g,d,info] = tvwcg1IQC(G,IQC,Tspan,NormType,Opt)
% NormType = 'L2toL2' or 'L2toE'. This input is optional with
% default 'L2toL2'.


%% Input Processing

% Check # of inputs
narginchk(3,5);
nin = nargin;
nout = nargout;

% Parse inputs
NormType = [];
Opt = [];
if nin==4
    if isa(varargin{end},'tvnormOptions')
        Opt = varargin{1};
    else
        NormType = varargin{1};
    end
elseif nin==5
    [NormType,Opt] = deal(varargin{:});
end

% Default Inputs
if isempty(NormType)
    NormType = 'L2toL2';
end
if isempty(Opt)
    Opt = tvnormOptions;
end

% Process Time Span
if isscalar(Tspan)
    T0=0;
    Tf=Tspan;
else
    T0 = Tspan(1);
    Tf = Tspan(2);
end

% Process Options
if isequal(Opt.Display,'on')
    DispFlag = true;
else
    DispFlag = false;
end
OdeOpt = tvodeOptions('OdeSolver',Opt.OdeSolver,'OdeOptions',Opt.OdeOptions);

%% IQC Information and Signal Dimensions
% XXX - Add error checking, e.g. if Delta is Nw-by-Nv then Psi is
% Nz-by-(Nv+Nw), M=M' is Nz-by-Nz, and blk is [Nw Nv].
% XXX - For now assume Psi is LTI and M is constant. These could be TV.

% Extract IQC
Psi = ss(IQC{1});
M = IQC{2};
blk = IQC{3};

% Signal Dimensions
Nv = blk(1);
Nw = blk(2);
Ne = size(G,1)-Nv;
Nd = size(G,2)-Nw;
%Nz = size(M,1);

%% Build Extended System

% State matrices for nominal system
[Ag,Bg,Cg,Dg] = ssdata(G);
Nx = size(Ag,1);
Bg1 = Bg(:,1:Nw);
Bg2 = Bg(:,Nw+1:end);
Cg1 = Cg(1:Nv,:);
Cg2 = Cg(Nv+1:end,:);

Dg11 = Dg(1:Nv, 1:Nw);
Dg12 = Dg(1:Nv, Nw+1:end);
Dg21 = Dg(Nv+1:end, 1:Nw);
Dg22 = Dg(Nv+1:end, Nw+1:end);

% State Matrices for IQC filter Psi
[Apsi,Bpsi,Cpsi,Dpsi] = ssdata(Psi);
Np = size(Apsi,1);
Bpsi1 = Bpsi(:,1:Nv);
Bpsi2 = Bpsi(:,Nv+1:end);
Dpsi1 = Dpsi(:,1:Nv);
Dpsi2 = Dpsi(:,Nv+1:end);

% State matrices for extended system
A = [Ag zeros(Nx,Np); Bpsi1*Cg1 Apsi];
B = [Bg1 Bg2; Bpsi1*Dg11+Bpsi2 Bpsi1*Dg12];
Ce1 = [Dpsi1*Cg1 Cpsi];
Ce2 = [Cg2 zeros(Ne,Np)];
De1 = [Dpsi1*Dg11+Dpsi2 Dpsi1*Dg12];
De2 = [Dg21 Dg22];

%% Build Cost Function Matrices
% Note:  R(t,g) = R0(t) - g^2 R1(t)
% E = [];
% Q = Ce1'*M*Ce1;
% R0 = De1'*M*De1;
% S = Ce1'*M*De1;
% R1 = blkdiag( zeros(Nw), eye(Nd) );
% if isequal(NormType,'L2toL2')
%     Q = Q + Ce2'*Ce2;
%     S = S + Ce2'*De2;
%     R0 = R0 + De2'*De2;
%     F = zeros(Nx+Np);
% else
%     % XXX Verify De2(t)==0 for L2 to E
%     % XXX This set-up requires R(t) to be inverted at each time step.
%     % However, no inversion is actually required since we can set
%     % R=I and B'=B/gamma.  This should speed up the integration.
%     Ce2Tf = tvsubs(Ce2,Tf);
%     F = Ce2Tf'*Ce2Tf;
% end

% State matrices
Cost.A = A;
Cost.B = B;

% IQC Cost
Qi =Ce1'*M*Ce1;
Ri = De1'*M*De1;
Si = Ce1'*M*De1;
Cost.IQC = struct('Q',Qi,'R',Ri,'S',Si);

% Performance Cost
R1p = blkdiag( zeros(Nw), eye(Nd) );
if isequal(NormType,'L2toL2')
    Qp = Ce2'*Ce2;
    Sp = Ce2'*De2;
    R0p = De2'*De2;
    Fp = zeros(Nx+Np);
else
    % XXX Verify De2(t)==0 for L2 to E
    % XXX This set-up requires R(t) to be inverted at each time step.
    % However, no inversion is actually required since we can set
    % R=I and B'=B/gamma.  This should speed up the integration.
    Qp = zeros(Nx+Np);
    Sp = zeros(Nx+Np,Nw+Nd);
    R0p = zeros(Nw+Nd);
    Ce2Tf = tvsubs(Ce2,Tf);
    Fp = Ce2Tf'*Ce2Tf;
end
Cost.Performance = struct('Q',Qp,'R0',R0p,'R1',R1p,'S',Sp,'F',Fp);

% Combined Cost
E = [];
Q = Qi + Qp;
S = Si + Sp;
R0 = Ri + R0p;
R1 = R1p;
F = Fp;
Cost.Combined = struct('Q',Q,'R0',R0,'R1',R1,'S',S,'F',F);


%% Lower Bound Phase
%  Require R(t,g)=R0(t)-g^2*R1(t)<0 for all t in [0,T]
% XXX This assumes the time grid is sufficiently fine.
if isequal(NormType,'L2toL2')
    ev = eig(R0,R1);
    ev = ev.Data(:);
    evmax = max( ev( isfinite(ev) ) );
    gLow = sqrt(evmax);
else
    gLow = 0;
end
gLow = max(gLow,Opt.Bounds(1));
if DispFlag
    fprintf('\n Lower Bound = %4.3f',gLow);
end

%% Upper Bound Phase
PUpp = []; PdotUpp = []; solUpp = [];
PLow = []; PdotLow = []; solLow = [];
if isfinite(Opt.Bounds(2))
    % User specified a (finite) upper bound.
    % Verify (or disprove) this upper bound.
    gTry = Opt.Bounds(2);
    R = R0-gTry^2*R1;
    if nout==1
        P = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
        Pdot = []; sol = [];
    else
        [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
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
        R = R0-gTry^2*R1;
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
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
AbsTol = Opt.AbsTol;
RelTol = Opt.RelTol;
if isfinite(gUpp)
    if DispFlag
        fprintf('\n');
    end
    while (gUpp - gLow > RelTol*gUpp+AbsTol)
        gTry = (gUpp+gLow)/2;
        
        % Solve LTV Riccati Equation
        R = R0-gTry^2*R1;
        if nout==1
            P = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
            Pdot = []; sol = [];
        else
            [P,~,Pdot,sol] = cdre(A,B,Q,R,S,E,F,[Tf T0],OdeOpt);
        end
        
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
g = [gUpp, gLow];

info.Lower.Gain = gLow;
info.Lower.P = PLow;
info.Lower.Pdot = PdotLow;
info.Lower.sol = solLow;

info.Upper.Gain = gUpp;
info.Upper.P = PUpp;
info.Upper.Pdot = PdotUpp;
info.Upper.sol = solUpp;

info.Cost = Cost;

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
    
    % Evaluate State/Cost Matrices on Same Time Grid as P
    R = tvmat( R0-gLow^2*R1 );
    [A,B] = evalt(A,B,T);
    if isequal(NormType,'L2toL2')
        [R,S] = evalt(R,tvmat(S),T);
    else
        S = tvmat(zeros(Nx,Nw+Nd));
    end
    
    % Simulate Transformed Hamiltonian system
    M = R\(PLow*B+S)';
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
    dData = cat(3,zeros(Nw+Nd,1,2), dData);
    d = tvmat(dData,dTime);
end