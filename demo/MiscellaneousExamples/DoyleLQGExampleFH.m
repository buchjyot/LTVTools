%% Doyle's Example - Finite Horizon Robustness Approach
%
% Reference: Doyle, J. (1978). Guaranteed margins for LQG regulators. IEEE
% Transactions on Automatic Control, 23(4), 756-757.

%% LTV systems parameters
%
%  dx/dt = A(t) x(t) + Bu(t) u(t) + Bp(t) p(t) + Bn(t) v(t)
%      y = C(t) x(t) + Du(t) u(t) + Dp(t) p(t) + Dn(t) v(t)

A = [1 1; 0 1];

Bp = [1; 1];
Bn = [0; 0];
Bu = [0; 1];
BBB = [Bp, Bn, Bu];

C = [1 0];

Dp = 0;
Dn = 1;
Du = 0;
DDD = [Dp, Dn, Du];

% Define Time Grid
T0 = 0;
Tf = 2;
Ts = 0.1;
TSpan = [T0,Tf];
TGrid = T0:Ts:Tf;

% Original System
P = ss(A,Bu,C,Du);
tvP = evalt(tvss(P),TGrid);

%% LQR-LQG Cost Data
% Continuous Time cost is E[ int_0^inf x'Qx + u'Ru dt ]
% Process and sensor noise are white noise, zero mean, Gaussian
% with variances W = E[w^2] and V = E[v^2].
q = 1000;
Q = q*[1;1]*[1 1];
R = 1;

sig = q;
Wcov = sig;
Vcov = 1;

% Options
OdeSolver = 'ode23s';
Bounds = [0,1e3];
RelTol = 1e-3;
AbsTol = 1e-4;
OdeOptions = odeset('RelTol',RelTol,'AbsTol',AbsTol);
tvopt = tvodeOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions);
tvnopt = tvnormOptions('OdeSolver',OdeSolver,'Bounds',Bounds,...
    'AbsTol',AbsTol,'RelTol',RelTol,'OdeOptions',OdeOptions);

%% LQR Design
% Assumes full-state available for feedback
Plqr = ss(A,Bu,eye(2),Du);
tvPlqr = evalt(tvss(Plqr),TGrid);

% LQR Gain: u(t) = -KlqrFH(t)*x(t)
[KlqrFH,PlqrFH] = tvlqr(tvPlqr,Q,R,[],[],TSpan,tvopt);

% Reshape the data and plot gains
figure;tvplot(KlqrFH,'LineWidth',2);title('LQR Gains');
xlabel('Time (s)');legend('K1','K2');grid on;

% Compare with analytical Infinite Horizon solution from Doyle
% Solution will only match if horizon is huge e.g. T = 1000
f = 2+sqrt(4+q);
Ka = [1 1]*f;
norm(KlqrFH.Data(:,:,1)-Ka)

%% Kalman Filter Design
% Kalman Filter System (Having n(t) as a Measurement Noise)
Pkf = ss(A,[Bu Bp],C,[Du Dp]);
tvPkf = evalt(tvss(Pkf),TGrid);

% Convert constant matrices to finite horizon tvmat
[A,Bu,C,Du] = ssdata(tvP);

% Kalman Filter Gain L
P0 = eye(2);
[~,LkfFH,PkfFH] = tvkalman(tvPkf,Wcov,Vcov,P0,TSpan);

% Reshape the data and plot gains
figure;tvplot(LkfFH,'LineWidth',2);title('Kalman Filter Gains');
xlabel('Time (s)');legend('L1','L2');grid on;

% Compare with analytical solution from Doyle
d = 2+sqrt(4+sig);
La = [1; 1]*d;
norm(LkfFH.Data(:,:,end)-La)

%% LQG Design
% Make the time vector consistent
KlqrFH = evalt(KlqrFH,TGrid);
LkfFH = evalt(LkfFH,TGrid);

% LQR Controller
Tlqr = feedback(tvPlqr,KlqrFH);

% Form LQG Controller
% Input is y and output is u
KlqgFH = tvss(A-Bu*KlqrFH-LkfFH*C,LkfFH,-KlqrFH,0);

% Form closed-loop
% Note: This is in positive feedback because the
% LQR gain already has the negative sign.
Tlqg = feedback(tvP*KlqgFH,1,+1);

% Compare closed-loop state-matrix with analytical expression from Doyle
Aa = [1 1 0 0; 0 1 -f -f; d 0 1-d 1; d 0 -d-f 1-f];
Acl = ssdata(Tlqg);
% norm(Aa-Acl.Data(:,:,end-500))

%% Uncertainty Weights
DelNorm = 1;
WDelSclin = ss(sqrt(DelNorm));
WDelSclout = ss(sqrt(DelNorm));

%% Disturbance Weights
Wd = ss(1);
Wn = ss(1);

%% LQR Uncertain Interconnection
tvPlqr = tvss(tvPlqr.A,[[1;1] tvPlqr.B],tvPlqr.C,[[0;0] tvPlqr.D]);
systemnames = 'tvPlqr KlqrFH WDelSclin WDelSclout Wd'; %#ok<*NASGU>
inputvar = '[w; d; p]';
outputvar = '[WDelSclin; tvPlqr(1)]';
input_to_tvPlqr = '[p; Wd+WDelSclout-KlqrFH]';
input_to_KlqrFH = '[tvPlqr]';
input_to_WDelSclin = '[Wd-KlqrFH]';
input_to_WDelSclout = '[w]';
input_to_Wd = '[d]';
cleanupsysic = 'yes';
Tlqrunc = sysic;

%% LQG Uncertain Interconnection
tvPf = evalt(tvss(A,BBB,[C;C],[zeros(1,3);DDD]),TGrid);
systemnames = 'tvPf KlqgFH WDelSclin WDelSclout Wd Wn';
inputvar = '[w; d; p; n]';
outputvar = '[WDelSclin;tvPf(1)]';
input_to_tvPf = '[p;Wn;Wd-KlqgFH+WDelSclout]';
input_to_KlqgFH = '[tvPf(2)]';
input_to_WDelSclin = '[Wd-KlqgFH]';
input_to_WDelSclout = '[w]';
input_to_Wd = '[d]';
input_to_Wn = '[n]';
cleanupsysic = 'yes';
Tlqgunc = sysic;

%% Nominal Closed-Loop
Tlqrnom = lft(0,Tlqrunc);
Tlqgnom = lft(0,Tlqgunc);

%% Save workspace to MAT file
save(mfilename);