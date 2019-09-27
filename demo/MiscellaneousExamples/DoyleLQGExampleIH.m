%% DoyleLQGExample
%
% Reference: Doyle, J. (1978). Guaranteed margins for LQG regulators. IEEE
% Transactions on Automatic Control, 23(4), 756-757.
%
% Paper showing that LQC regulators can have poor margins.

%% Plant Data
%  dx/dt = A x + Bu u + Bw w + Bv v
%      y = C x + Du u + Dw w + Dv v; 
A = [1 1; 0 1];

Bu = [0; 1];
Bw = [1; 1];
Bv = [0; 0];

C = [1 0];

Du = 0;
Dw = 0;
Dv = 1;

%% LQG Cost Data
% Cost is E[ int_0^inf x'Qx + u'Ru dt ]
% and process/sensor noise is Qw = E[w^2] and Qv = E[v^2];
q = 1000;
Q = q*[1;1]*[1 1];
R = 1;

sig = q; 
Qw = sig;   
Qv = 1;   

%% Solve LQG Problem 
% Use lqr/kalman to solve for the state feedback and estimator gains

% LQR Gain: u = -Klqr*x
[Klqr,Xlqr] = lqr( ss(A,Bu,eye(2),0), Q,R);

% Compare with analytical solution from Doyle
f = 2+sqrt(4+q);
Ka = [1 1]*f;
norm(Klqr-Ka)

% Kalman Filter Gain
[~,Lkf,Xkf] = kalman( ss(A,[Bu Bw],C,[Du Dw]),Qw,Qv);

% Compare with analytical solution from Doyle
d = 2+sqrt(4+sig);
La = [1; 1]*d;
norm(Lkf-La)

%% Compute Margins

% Form plant (without noises)
% Input is u and output is y
G = ss(A,Bu,C,0);

% Form LQG Controller
% Input is y and output is u
Clqg = ss(A-Bu*Klqr-Lkf*C,Lkf,-Klqr,0);

% Form closed-loop 
% Note: This is in positive feedback because the
% LQR gain already has the negative feedback.
CL = feedback(G*Clqg,1,+1);

% Compare closed-loop state-matrix with analytical expression from Doyle
Aa = [1 1 0 0; 0 1 -f -f; d 0 1-d 1; d 0 -d-f 1-f];
norm(CL.a-Aa)

% Classical and disk margins
[cm,dm] = loopmargin(G,-Clqg);
cm.GainMargin
dm

% Bode plots of plant G, LQG controller -Clqg, and loop L
% Note: Include negative sign in Clqg to correspond to more standard
% negative feedback loop.
L = G*-Clqg;

figure(1)
bode(G,'b',-Clqg,'r--',L,'g--');

% Limiting controller as q=sig-->inf (?)
Clim = tf([-2 1],1);
Llim = G*-Clim;
hold on;
bode(-Clim,'c',Llim,'c');
hold off

return
%% Solve LQG Problem 
% Use H2SYN to compute the optimal controller. This is simply
% for comparison. The construction of the "error" assumes that w 
% and v are scalar signals.

Bh2 = [Bw*sqrt(Qw) Bv*sqrt(Qv) Bu];
Ch2 = [sqrt(q)*[1 1]; 0 0; C];
D = [0 0 0; 0 0 sqrt(R); Dw*sqrt(Qw) Dv*sqrt(Qv) Du];
P = ss(A,Bh2,Ch2,D);
[Kh2,CLh2,GAMh2,INFOh2] = h2syn(P,1,1);

figure(2);
bode(Clqg,'b',Kh2,'r--')
