%% Horizon
T0 = 0;
Tf = 3;

%% Systems
plant_case = 1;
switch plant_case
    case 1
        % SISO Plant
        Nx = 2;
        Ny = 1;
        Nu = 1;
        G = ss(tf(1,[1 0.6 1]));
        
    case 2
        % SIMO Plant
        Nx = 2;
        Ny = 1;
        Nu = 1;
        G = ss(tf(1,[1 0.6 1]));
        G = ss(G.A,G.B,eye(2),0);
        
    case 3
        % MIMO Plant
        rng(0);
        Nx = 4;
        Nu = 2;
        Ny = 3;
        G = rss(Nx,Ny,Nu);
        G.D = 0;
        
    case 4
        % Example in which the feedthrough matrix has two identical
        % singular values i.e. there are two input directions for which you
        % get the same amplification. In this case, the computed
        % disturbances using RDE will not match up with power iterations.
        rng(0);
        D = randn(3);
        [U,S,V] = svd(D);
        S(2,2) = S(1,1);
        D = U*S*V';
        zeta = 0.2;
        wn = 4;
        G = tf([2*zeta*wn  0],[1 2*zeta*wn wn^2])*D;
        [Ny,Nu] = size(G);
        
    case 5
        % Static Gain
        G = ss(rand(5));
        [Ny,Nu] = size(G);
end

%% Power Iterations
NE = 0;
Gt = tvss(G,[T0,Tf]);
pOpt = poweritOptions('Display','on','StoreAllIter',true);
pSpec= poweritSignalSpec('NE',NE);
[g1,d1,info] = powerit(Gt,[T0,Tf],pSpec,pOpt);
% analyzeResults(info)

%% RDE Bisections
tOpt = tvnormOptions('Display','on');
[g2,d2] = tvnorm(Gt,tOpt);

%% Plot Disturbances
figure(1);
tvplot(d1,'b',d2,'--r','LineWidth',2);
grid on; box on;
ylabel('Disturbances');