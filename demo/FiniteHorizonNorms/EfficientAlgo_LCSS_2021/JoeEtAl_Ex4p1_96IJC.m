%% LTV System
% This example uses the Example 4.1 in 
% 
% Imae, J. and Wanyoike, G., 1996. Hinf norm computation for LTV systems
% using nonlinear optimal control algorithms. International Journal of
% Control, 63(1), pp.161-182.

t = 0:0.1:10;
At = [tvmat(-1 + sin(t),t) 1; 0 -4]; 
Bt = eye(2);
Ct = eye(2);
G = tvss(At,Bt,Ct,0);

% Options
pOpt   = poweritOptions('StopTol',1e-2/5,'Display','off');
tvnOpt = tvnormOptions('AbsTol',1e-2,'Display','off');
tvOpt  = tvodeOptions;

% Initial Disturbance
rng(0);
Nt = length(t);
d0 = tvmat(randn(2,1,Nt),t);
d0 = d0/tvnorm(d0);

%% Power Method
[gPower,dPower,infoPower] = tvnormPower(G,0,d0,[],pOpt);

%% Bisection
[gRDE,dRDE,infoRDE] = tvnormBisect(G,0,tvnOpt);

%% Combined Method
[gC,dC,infoC] = tvnormCombined(G,0,d0,tvnOpt);
save(mfilename);
return;

%% Plot disturbances
figure(1);clf;
tvplot(-dC,'LineWidth',2);
legend('d_1(t)','d_2(t)','Location','northwest');grid on;box on;
ylabel('Worst-Case Disturbances');
xlabel('Time (sec)')