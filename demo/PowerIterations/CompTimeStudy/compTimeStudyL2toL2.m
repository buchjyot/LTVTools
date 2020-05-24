%% Power Iteration Example
% This file compares the computational time as state dimention is varied
% for computing finite horizon worst-case disturbance and induced L2 to L2
% gain lower bound using the following two methods.
%
%   a) Power iteration
%   b) Riccati approach usd in tvnorm function

%% Create MAT file for the Plant Data
if false
    % We use the following state dimentions for this study.
    NxAll = [1 10:10:150];
    m = length(NxAll);
    nPlant = 5;
    
    % Set seed
    rng(0);
    
    %  We create 5 random SISO plants using rss for a given state dimention and
    % average the computational time to represent a sample.
    k = 0;
    for i = 1:m
        for j = 1:nPlant
            k = k + 1;
            Gall(:,:,k) = rss(NxAll(i),1,1);
        end
    end
    
    % Save MAT file
    filename = sprintf('LTIDataSet_nPlant%d',nPlant);
    save(filename,'Gall','k','nPlant','NxAll');return;
else
    % Load MAT file
    load('LTIDataSet_nPlant5');
end

%% Fixed Horizon
T0 = 0;
Tf = 3;

%% Memory Allocation
gLB     = zeros(k,1);
tvnLB   = zeros(k,1);
tP      = zeros(k,1);
tB      = zeros(k,1);
nBisect = zeros(k,1);
nPiter  = zeros(k,1);

%% Options
tvpOpt = tvpoweritOptions('StopTol',1e-2);
tvnOpt = tvnormOptions('RelTol',1e-2,'AbsTol',1e-2);

%% Main For Loop
for i = 1:k
    % Time-varying representation
    G = evalt(tvss(Gall(:,:,i)),[T0,Tf]);
    
    % Power iterations
    [gLB(i),~,info1] = tvpowerit(G,tvpOpt);
    tP(i) = info1.TotalTime;
    nPiter(i) = info1.TotalIter;
    
    % TVNORM
    [tvn,~,info2] = tvnorm(G,tvnOpt);
    tvnLB(i) = tvn(1);
    nBisect(i) = info2.TotalBisections;
    tB(i) = info2.TotalTime;
    
    % Display iteration number
    disp(i);
end

% Save Data
filename = [mfilename sprintf('_nPlant%d',nPlant)];
save(filename,'gLB','nBisect','nPiter','tB','tP','tvnLB','nPlant','NxAll');

%% Plot Data
load('compTimeStudyL2toL2_nPlant5.mat');

% Average Per Model Order
k = 1;
for i = 1:nPlant:length(tP)
    nBisect_AVG(k) = sum(nBisect(i:i+nPlant-1))/nPlant; %#ok<*SAGROW>
    nPiter_AVG(k)  = sum(nPiter(i:i+nPlant-1))/nPlant;
    tB_AVG(k)      = sum(tB(i:i+nPlant-1))/nPlant;
    tP_AVG(k)      = sum(tP(i:i+nPlant-1))/nPlant;
    k = k + 1;
end

% Plot Computational Time
figure(1);clf;
plot(NxAll,tB_AVG,'--bo',NxAll,tP_AVG,'--rs','LineWidth',2);
grid on;box on;
xlabel('System Order (Nx)');
ylabel('Average Computational Time (sec)');
legend('Bisection Method','Power Iteration','Location','northwest');

% Plot Average Number of Iterations
figure(2);clf;
plot(NxAll,nBisect_AVG,'--bo',NxAll,nPiter_AVG,'--rs','LineWidth',2);
grid on;box on;
xlabel('System Order (Nx)');
ylabel('Average Number of Iterations');
legend('Bisections','Power Iterations');

% Maximum error in lower bound
max(abs((gLB - tvnLB)./max(gLB,tvnLB)))