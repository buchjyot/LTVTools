%% System
% Horizon
T0 = 0;
Tf = 3;
G = ss(tf(1,[1 0.6 1]));
C = eye(2);
G = ss(G.A,[],C,[]);
Gt = tvss(G,[T0,Tf]);

% Initial condition ellipsoid
% rng(2357);
% M = randn(2);
% E0 = M'*M;
E0 = [0.1099, -0.1312; -0.1312, 0.9807];
E0 = E0/norm(E0);

% Options
pOpt = poweritOptions('Display','on','StoreAllIter',true);
pSpec= poweritSignalSpec('NE',2,'InitialConditions','free','InitialCondCostMat',E0);

% power iterations
[g,dwc,info] = powerit(Gt,[T0,Tf],pSpec,pOpt);

%% Plot
figure(1);clf;box on;hold on;grid on;

Xwc = info.Xwc;
Xwc0 = tvsubs(Xwc,T0);
XwcT = tvsubs(Xwc,Tf);
plot(Xwc(1),Xwc(2),'r','LineWidth',2);
plot(Xwc0(1),Xwc0(2),'ko','MarkerFaceColor','k');
plot(XwcT(1),XwcT(2),'ro','MarkerFaceColor','r');

p1 = plot_gaussian_ellipsoid([0,0],inv(E0));
patch(p1.XData,p1.YData,'b');
alpha(0.1);

p2 = plot_gaussian_ellipsoid([0,0],g^2*(C'*C));
patch(p2.XData,p2.YData,'r');
alpha(0.1);

xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);axis equal