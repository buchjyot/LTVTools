%% Load Data
load('MSDEx_RobustAnalysis');
Gt = tvss(Gunc,[0,Tall(1)]);

%% Verify worst-case gain for computed input
Uwc = wcSig{1};

% wcgSim should match with wcgLB
[wcgSim,Ysim,Xsim] = tvlsimwcgain(Gt,Uwc,1,1,NE);

%% Read worst-case signals and check causality
% Read inputs
w = Uwc(1);
u = Uwc(2);
Usim = Uwc;

% Read worstcase outputs
v = Ysim(1);
yL2 = Ysim(2:3);

% Compute norms
tvnw = tvnorm(w);
tvnu = tvnorm(u);
tvnv = tvnorm(v);

% Plot signals
figure(1)
tvsplot(Ysim,'b','LineWidth',2);
subplot(3,1,1);ylabel('v');grid on; box on;
subplot(3,1,2);ylabel('e_1');grid on; box on;
subplot(3,1,3);ylabel('e_2');grid on; box on;

figure(2),tvsplot(Usim,'r','LineWidth',2);
subplot(2,1,1);ylabel('w');grid on; box on;
subplot(2,1,2);ylabel('d');grid on; box on;

% Unify time grid
tUnion = union(w.Time,v.Time);
[v,w] = evalt(v,w,tUnion);

% Integrate state dynamics for delta
fvw = v'*v - w'*w;
x_delta_dot = @(t,x) tvsubs(fvw,t);
odeopt = odeset('RelTol',1e-6,'AbsTol',1e-6);
[tode,x_delta] = ode23s(x_delta_dot,[0 Tall(1)],0,odeopt);
x_del = tvmat(x_delta,tode);

% Plot States
figure(3);clf;
tvplot(x_del,'b','LineWidth',2);
ylabel('x_\Delta(t)');grid on;box on;

% x_del must be negative to meet causality requirement

%% Construct Specific Bad Nonlinearity
% Consider w = Delta(v)
DeltaBad = tvnw/tvnv; % wcuIH.Del

% Take the plant and wrap uncertainty
Fu = lft(DeltaBad,Gt);

% Compute the noninal input-output performance
tvOpt = tvnormOptions('Display','on');
[gbnd,dwc] = tvnorm1(Fu,NE,tvOpt);