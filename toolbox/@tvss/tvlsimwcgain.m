function [wcg,Y,X] = tvlsimwcgain(G,U,Nv,Nw,NE,x0,Opt)
%% TVLSIMWCGAIN
%
% Compute induced gain for the uncertain plant given specific bad input U

% Make sure the horizon is real number
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);
Nx = order(G);
[NY,NU] = size(G);
Nu = NU - Nw;

% Input Processing
narginchk(4,7);
nin = nargin;
switch nin
    case 4
        NE = 0;
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
    case 5
        x0 = zeros(Nx,1);
        Opt = tvodeOptions;
    case 6
        Opt = tvodeOptions;
end

% Check if NE is valid
if ~(NE >= 0) || ~(NE <= NY-Nv)
    error('Euclidean penalty must be between 0 and total number of plant outputs NY-Nv.');
end
NL2 = NY-Nv-NE;

% Input Signals
w = U(1:Nw,:);
u = U(Nw+1:Nw+Nu,:);
U1 = [w;zeros(Nu,1)];
U2 = [zeros(Nw,1);u];

% Linear simulation
[Y1,X1] = tvlsim(G,U1,[T0,Tf],x0,Opt);
[Y2,X2] = tvlsim(G,U2,[T0,Tf],x0,Opt);

% Bisect on sqrt(alpha)
% Find alpha (or 1/wcg) such that tvnorm of output is 1

% Evaluate on same time grid
tUnion = union(Y1.Time,Y2.Time);
[Y1,Y2,X1,X2] = evalt(Y1,Y2,X1,X2,tUnion);

% Construct a1,a2,b1,b2,c1,c2 such that
% v   = a1 + a2*sqrt(alpha)
% yL2 = sqrt(alpha)*(b1 + b2*sqrt(alpha))
% yE  = sqrt(alpha)*(c1(T) + c2(T)*sqrt(alpha));

a1 = Y1(1:Nv,:);
b1 = Y1(Nv+1:Nv+NL2,:);
c1 = tvsubs(Y1(Nv+NL2+1:Nv+NL2+NE,:),Tf);

a2 = Y2(1:Nv,:);
b2 = Y2(Nv+1:Nv+NL2,:);
c2 = tvsubs(Y2(Nv+NL2+1:Nv+NL2+NE,:),Tf);

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
        [sa,wcg] = LOCALBisect(a1,a2,b1,b2,c1,c2,gLowIn,gUppIn);
    end
else
    [sa,wcg] = LOCALBisect(a1,a2,b1,b2,c1,c2,0,1e6);
end

% Due to linearity perform super-position
X = X1 + X2*sa;
OutScl = blkdiag(eye(Nv),sa*eye(NL2+NE));
Y = OutScl*(Y1 + Y2*sa);
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