%% Plant Data
% Horizon
T0 = 0;
Tf = 3;

% Problem data
Nx = 2;
Ny = 1;
Nu = 1;
G = ss(tf(1,[1 0.6 1]));
G = ss(G.A,G.B,eye(2),0);
F = G.C'*G.C;

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

% Options
NE = 2;
tOpt = tvnormOptions('Display','on','Bounds',[0 10]);

% Uncertainty in Initial Condition
E0_CASE = 4;
switch E0_CASE
    case 1
        rng(7582);
        Er = rand(2);
        E0 = Er'*Er;
        E0 = E0/norm(E0);
    case 2
        E0 = eye(NE);
    case 3
        E0 = diag(inf(2,1));
    case 4
        E0 = [0.1099, -0.1312; -0.1312, 0.9807];
end

%% RDE Approach
[g1,d1,info1] = tvnormic(Gt,E0,NE,tOpt);
P = info1.Upper.P;
P0 = tvsubs(P,T0);
nd1 = tvnorm(d1);

%% Power Iteration Lower Bound
pOpt  = poweritOptions('Display','on','StopTol',1e-2,'OdeSolver','ode23s','StoreAllIter',true);
pSpec = poweritSignalSpec('NE',NE,'InitialConditions','free','InitialCondCostMat',E0);

% Power Iterations
[g2,d2,info2] = powerit(Gt,[T0,Tf],pSpec,pOpt);
analyzeInfo(info2);
nd2 = tvnorm(d2);

%% Saperate figure for P0 and E0
figure(5);clf;

% Plot P(0) Ellipsoid
p1 = plot_gaussian_ellipsoid([0,0],P0);
p1.Color = 'b';grid on;box on;hold on;
pp1 = patch(p1.XData,p1.YData,'b'); alpha(0.1);

% Plot E0 Ellipsoid
p2 = plot_gaussian_ellipsoid([0,0],E0*g1(1)^2);
p2.Color = 'g';
pp2 = patch(p2.XData,p2.YData,'g'); alpha(0.1);

% Finalize plot
axis equal
% legend([p1 p2],{'P(0)','\gamma_{ub}^2 E0'});
xlabel('x_1','FontSize',12);
ylabel('x_2','FontSize',12);

%% Plot Ellipsoids

% Ellipsoids
figure(6);clf;grid on;box on;hold on;

% Plot inv(E0) Ellipsoid
p3 = plot_gaussian_ellipsoid([0,0],inv(E0));
p3.Color = 'b';
pp3 = patch(p3.XData,p3.YData,'b'); alpha(0.1); %#ok<*NASGU>

% Plot trajectory in state-space
Xwc = info2.Xwc;
X0  = tvsubs(Xwc,T0);
tvplot(Xwc(1),Xwc(2),'r','LineWidth',2);
plot(X0(1),X0(2),'ko','LineWidth',2,'MarkerFaceColor','k');

% % Plot a vector pointing to the initial conditions
% quiver(0,0,X0(1),X0(2),1,'m','LineWidth',1.5);
% plot(0,0,'mo','LineWidth',1.5,'MarkerFaceColor','m');

% Highlight final conditions
XT = tvsubs(Xwc,Tf);
plot(XT(1),XT(2),'ro','LineWidth',2,'MarkerFaceColor','r');

% Plot final state ball
p4 = plot_gaussian_ellipsoid([0,0],F*g2^2);
p4.Color = 'r';
pp4 = patch(p4.XData,p4.YData,'r'); alpha(0.1);

% Finalize plot
axis equal
%legend([p3 p4],{'E0^{-1}','\gamma_{lb}^2 '});
xlabel('x_1','FontSize',12);
ylabel('x_2','FontSize',12);

%% Plot worst-case disturbances
% Normalize Disturbances
d1 = d1/nd1;
d2 = d2/nd2;

figure(7);clf;
tvplot(d1,'b',d2,'r--','LineWidth',2);
title('Worst-case disturbances');grid on;box on;
legend('d_{RDE}','d_{POW}');