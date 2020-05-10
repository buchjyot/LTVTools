%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

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
        % MIMO Plant
        rng(0);
        Nx = 20;
        Nu = 2;
        Ny = 3;
        G = rss(Nx,Ny,Nu);
end

% Time-varying representation
G   = evalt(tvss(G),[T0,Tf]);

%% Power Itertaions
% Perform iterations
fprintf('### Power Iterations:\n');
pOpt = poweritOptions('Display','on');
[gLB,dwc1,info1] = powerit(G,pOpt);

% Display
fprintf('L2 Gain Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB,info1.TotalIter,info1.TotalTime);

%% Verify with LTVTools Lower Bound
[tvn,dwc2,info2] = tvnorm(G);

% Display
fprintf('### Riccati Approach :\n');
fprintf('L2 Gain Lower Bound: %.3f, Total RDE Bisections: %d, Computational Time: %.3f seconds\n',tvn(1),info2.TotalBisections,info2.TotalTime);

%% Plot worst-case disturbances
figure(1);clf;hold on;
fh1 = tvplot(dwc1,'b','LineWidth',2.5);
fh2 = tvplot(-dwc2,'r--','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Wort-Case Disturbance Signal (dwc)');
legend('Power Iteration','Riccati Bisection');