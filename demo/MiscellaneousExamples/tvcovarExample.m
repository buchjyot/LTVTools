%% tvcovarExample
% This example tests the tvcovar function to compute the covariance of the
% time-varying state space systems.

%% LTV System
G = ss(tf(1,[1 0.2 1]));
SYS = ss(G.A,G.B,eye(2),[G.D;G.D]);
SYSt = tvss(SYS);

% Large Horizon
T0 = 0;
Ts = 0.1;
Tf = 50;

%% Run TVCOVAR
NoiseLevel = 1;
[Pt,Qt] = tvcovar(evalt(SYSt,T0:Ts:Tf),NoiseLevel);
[P,Q] = covar(SYS,NoiseLevel);

%% Plot
figure(1);
tvplot(Qt,'b',evalt(tvmat(Q),Qt.Time),'--r');
title('State Covariances');

figure(2);grid on;box on;
Tall = [0.3,0.5,1,2,3,10,20];
NT = length(Tall);
legendArray = cell(NT,1);
for i = 1:NT
    % Plot 1 sigma (standard deviation) ellipse for states
    h = plot_gaussian_ellipsoid([0,0],tvsubs(Qt,Tall(i)),1);
    h.LineWidth = 2;
    legendArray{i} = sprintf('Tf = %.1f',Tall(i));
end
plot(0,0,'ok','MarkerFaceColor','k')
legend(legendArray,'Location','bestoutside');
axis equal;
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
title('1-\sigma ellipses around initial condition [0,0]');