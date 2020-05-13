%% Horizon
T0 = 0;
Tf = 3;

%% Systems
plant_case = 2;
switch plant_case
    case 1
        % SISO Plant
        Nx = 2;
        Ny = 1;
        Nu = 1;
        G = ss(tf(1,[1 0.6 1]));
        
    case 2
        % MIMO Plant
        rng(0);
        Nx = 4;
        Nu = 2;
        Ny = 3;
        G = rss(Nx,Ny,Nu);
        
    case 3
        % Example in which the feedthrough matrix has two identical
        % singular values i.e. there are two input directions for which you
        % get the same amplification. Thus, power iteration method may get
        % stuck, because input directions keeps bouncing.
        rng(0);
        D = randn(3);
        [U,S,V] = svd(D);
        S(2,2) = S(1,1);
        D = U*S*V';
        zeta = 0.2;
        wn = 4;
        G = tf([2*zeta*wn  0],[1 2*zeta*wn wn^2])*D;
        [Ny,Nu] = size(G);
        
    case 4
        % Static Gain
        G = ss(rand(5));
        [Ny,Nu] = size(G);
end

%% Nominal SimIn
time = [T0;Tf];
Nt = length(time);
data = zeros(Nt,Nu);
simin = [time data];