%% Constant D-Scale Example

%% Uncertain System
% Uncertain system is Fu(Delta,G).  Assume Delta is a memoryless function
% in the sector [-1 1] for this example.
load LTIUncExData1;
G = MBad(:,:,7);
[Ag,Bg,Cg,Dg] = ssdata(G);

Ag = round(Ag*10)/10;
Bg = round(Bg*10)/10;
Cg = round(Cg*10)/10;
Dg = round(Dg*10)/10;
G = ss(Ag,Bg,Cg,Dg);

%% Infinite Horizon Performance vs Constant D-scale

% Use LMI condition to bound the infinite-horizon L2 gain
nx = size(Ag,1);
cvx_begin sdp
    % Create a symmetric 2-by-2 matrix variable
    variable P(nx,nx) symmetric;
    variable gsq;
    variable lam;
    
    % LMI #1 
    P >= 0; 
    
    % LMI #2
    Riqc = [Cg(1,:) Dg(1,:); zeros(1,nx) 1 0];
    Miqc = [1 0; 0 -1];
    
    Rg = [Cg(2,:) Dg(2,:); zeros(1,nx) 0 1];
    Mg = [1 0; 0 -gsq];
    [Ag'*P+P*Ag P*Bg; Bg'*P zeros(2)] + lam*Riqc'*Miqc*Riqc+Rg'*Mg*Rg <=0;
    
    % Objective function
    minimize( gsq );
cvx_end

% Print Results
lamopt = lam;
dopt = sqrt(lamopt);
gopt = sqrt(gsq);
fprintf('\n Gain bound = %4.4f \t Optimal D-scale = %4.4f\n',gopt,dopt);

Lfac = diag([dopt 1]);
Rfac = diag([dopt gopt]);

if false
    % Verify LMI
    LMI1 = [Ag'*P+P*Ag P*Bg; Bg'*P zeros(2)] + lam*Riqc'*Miqc*Riqc+Rg'*Mg*Rg;
    eig(LMI1)
    
    % Rewrite LMI
    fac1 = diag([sqrt(lam) 1])*[Cg Dg];
    fac2 = diag([-lam,-gsq]);
    LMI2 = [Ag'*P+P*Ag P*Bg; Bg'*P fac2] + fac1'*fac1;
    eig(LMI2)
    
    % Transform LMI to Bounded Real Form
    Bscl = Bg/Rfac;
    Cscl = Lfac*Cg;
    Dscl = Lfac*Dg/Rfac;
    LMI3 = [Ag'*P+P*Ag P*Bscl; Bscl'*P -eye(2)] + [Cscl Dscl]'*[Cscl Dscl];
    eig(LMI3)
end

% Verify norm of Scaled System
n = norm( Lfac*G/Rfac, inf);
fprintf(' Norm of Scaled System = %4.4f (Should be <1)\n',n);

% Graphically display norm of scaled system vs. D-scale
Nd = 250;
dscale = linspace(1,3,Nd);
Jih = zeros(Nd,1);
for i=1:Nd
    Dl = [dscale(i) 0; 0 1];
    Dr = [1/dscale(i) 0; 0 1];
    S = [1 0; 0 1/gopt];
    Jih(i) = norm( Dl*G*S*Dr, inf,1e-5);
end
    
figure(1); clf
plot(dscale,Jih,'b',dopt,1,'ro');
xlabel('D-scale');
ylabel('Scaled Gain');
[Jmin,Jidx] = min(Jih);
title( sprintf('min J = %4.4f at d= %4.4f ',Jmin,dscale(Jidx)) );
grid

%% Performance vs. Finite Horizon Time

% Create TVSS Object
Gtv = tvss(G);

% Time Horizon
Tdata = [10 15 20 30 40 60 90];
NT = numel(Tdata);

% IQC (Using optimal D-Scale)
Psi = ss( dopt*eye(2) );
M = [1 0;0 -1];
blk = [1 1];
IQC = {Psi,M,blk};

% Compute gain bound on various finite horizons
Jfh = zeros(NT,1);
for i=1:NT
    [g,d,info] = tvwcg1IQC(Gtv,IQC,Tdata(i));
    Jfh(i) = g(1);
    [i Tdata(i) Jfh(i)]
end

% Plot results
figure(2); clf
plot(Tdata,Jfh,'b',Tdata([1 end]),gopt*[1 1],'r--');
xlabel('Finite Horizon, T');
ylabel('Gain Bound');
xlim( Tdata([1 end]) )

%% Performance vs. D-scale (At One Finite Horizon Time)
T = 20;    % T=150;
Psi = ss( eye(2) );
M = dopt^2*[1 0;0 -1];
IQC = {Psi,M,blk};
[gT,dT,infoT] = tvwcg1IQC(Gtv,IQC,T);

if T==20
    dscale = unique(dopt*[1 linspace(0.89,1,5) 1.1 2]);
    %dscale = unique([dopt linspace(1.89,1.95,20) 2 2.1 2.3]);
else
    dscale = dopt*[linspace(0.99,1.02,10) 1.1 2 4];
end

Nd = numel(dscale);
Jdfh = zeros(Nd,2);
Jdih = zeros(Nd,1);
for i=1:Nd
    Psi = ss( eye(2) );
    M = dscale(i)^2*[1 0;0 -1];
    blk = [1 1];
    IQC = {Psi,M,blk};
        
    [g,d,info] = tvwcg1IQC(Gtv,IQC,T);
    Jdfh(i,:) = g;    
    
    cvx_begin sdp quiet
        % Create a symmetric 2-by-2 matrix variable
        variable P(nx,nx) symmetric;
        variable gsq;
    
        % LMI #1
        P >= 0;

        % LMI #2
        lam = dscale(i)^2;
        Riqc = [Cg(1,:) Dg(1,:); zeros(1,nx) 1 0];
        Miqc = [1 0; 0 -1];

        Rg = [Cg(2,:) Dg(2,:); zeros(1,nx) 0 1];
        Mg = [1 0; 0 -gsq];
        [Ag'*P+P*Ag P*Bg; Bg'*P zeros(2)] + lam*Riqc'*Miqc*Riqc+Rg'*Mg*Rg <=0;

        % Objective function
        minimize( gsq );
    cvx_end
    Jdih(i) = sqrt(gsq);
    
    [i dscale(i) Jdfh(i) Jdih(i)]
end
    
[Jmin,Jidx] = min(Jdfh(:,1));
idx = find( dscale==dopt );

figure(3); clf
plot(dscale,Jdfh,'b',dscale,Jdfh,'bx',...
    dscale,Jdih,'r--',dscale,Jdih,'rx',...
    dopt,gT,'co',dscale(idx),Jdih(idx),'ko');
xlabel('D-scale');
ylabel('Gain');
title( sprintf('min J = %4.4f at d= %4.4f ',Jmin,dscale(Jidx)) );
grid


%% Construct Cutting Plane
%dscale = 1.8; %1.9009; 
dscale= 3;
Psi = ss( eye(2) );
M = dscale^2*[1 0;0 -1];
IQC = {Psi,M,blk};
opt.RelTol = 0; opt.AbsTol = 2; opt.Display = 'on'; opt.Bounds = [4 7];
[g,d,info] = tvwcg1IQC(Gtv,IQC,T,opt);

A = info.Cost.A;
B = info.Cost.B;
[yy,xx] = tvlsim( evalt(tvss(A,B,eye(size(A)),zeros(size(B))),d.Time), d);
x = evalt(xx,d.Time);

Q = info.Cost.Combined.Q;
S = info.Cost.Combined.S;
R0 = info.Cost.Combined.R0;
R1 = info.Cost.Combined.R1;
R = R0 -  g(2)^2*R1;
J = [x;d]'*[Q S; S' R]*[x; d];

trapz(J.Time, J.Data(:))

Qi = info.Cost.IQC.Q; Ri = info.Cost.IQC.R; Si = info.Cost.IQC.S;
Qp = info.Cost.Performance.Q; Sp = info.Cost.Performance.S;
R0p = info.Cost.Performance.R0; R1p = info.Cost.Performance.R1;
Fp = info.Cost.Performance.F;

a0 = [x;d]'*[Qp Sp; Sp' R0p]*[x;d];
a1 = [x;d]'*[Qi Si; Si' Ri]*[x;d];
a2 = d'*R1p*d;

a0 = trapz(a0.Time,a0.Data(:));
a1 = trapz(a1.Time,a1.Data(:));
a2 = trapz(a2.Time,a2.Data(:));
tmp = a0+a1-g(2)^2*a2

dd = linspace(1.5,4.5,50);
gg = sqrt( (a0+a1*(dd/dscale).^2)/a2 );
hold on;
ph = plot(dd,gg,'c');
hold off