%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

%% Horizon
T0 = 0;
Tf = 5;

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
        % MIMO Plant
        rng(0);
        Nx = 20;
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
        
    case 4
        % Static Gain
        G = ss(rand(5));
end

% Time-varying representation
G   = evalt(tvss(G),[T0,Tf]);

%% Power Itertaions
% Perform iterations
fprintf('### Power Iterations:\n');
pOpt = poweritOptions('Display','on','StepSize','Auto');
[gLB,dwc1,info1] = powerit(G,pOpt);

% Display
fprintf('L2 Gain Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB,info1.TotalIter,info1.TotalTime);

%% Verify with LTVTools Lower Bound
[tvn,dwc2,info2] = tvnorm(G);

% Display
fprintf('### Riccati Approach :\n');
fprintf('L2 Gain Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n',tvn(1),info2.TotalBisections,info2.TotalTime);

%% Plot worst-case disturbances
figure(1);clf;hold on;box on;
fh1 = tvplot(dwc1,'b','LineWidth',2.5);
fh2 = tvplot(-dwc2,'r--','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1)],{'Power Iteration','Riccati Bisection'});