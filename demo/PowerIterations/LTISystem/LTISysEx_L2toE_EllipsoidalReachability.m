%% Plant Data
% Horizon
T0 = 0;
Tf = 3;

% Problem data
Nx = 2;
Ny = 1;
Nu = 1;
G = ss(tf(1,[1 0.6 1]));
C = eye(2);
G = ss(G.A,G.B,C,0);

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

%% Solve LDE
[A,B] = ssdata(Gt);
P = cdle(A,B,inv(E0),[T0,Tf]);
Po = C*P*C';
g = sqrt(max(eig(Po)));

figure(1);clf;
tvplot(g,'LineWidth',2);
grid on; box on;
ylabel('Euclidean Gain Bound');