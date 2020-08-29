%% Power Iteration Example
% This file compares only single iteration of the following methods.
%
%   a) Power iteration
%   b) Riccati approach usd in tvnorm function

%% Create MAT file for the Plant Data
% Load MAT file
load('LTIDataSet_nPlant5');

%% Fixed Horizon
T0 = 0;
Tf = 15;

%% Memory Allocation

% Power Method
pGain               = zeros(k,1);
pTotalTime          = zeros(k,1);
pTotalIter          = zeros(k,1);
pTimeAllIter        = cell(k,1);
pSinglePowerIt      = zeros(k,1);

% RDE Bisections
rUB                 = zeros(k,1);
rLB                 = zeros(k,1);
rTotalBisections    = zeros(k,1);
rTotalTime          = zeros(k,1);
rSingleRDE          = zeros(k,1);

% Combined method
cUB                 = zeros(k,1);
cLB                 = zeros(k,1);
cRDECount           = zeros(k,1);
cPowerIterations    = cell(k,1);
cTotalTime          = zeros(k,1);

%% Induced L2 gian
NE  = 0;
NL2 = 1;
Ny  = 1;
Nu  = 1;

%% Initial Disturbance
Tgrid = T0:1:Tf;
rng(0);
Nt = length(Tgrid);
d1 = tvmat(randn(Nu,1,Nt),Tgrid);
d1 = d1/tvnorm(d1);

%% Options
pOpt   = poweritOptions('StopTol',1e-2/5);
tvnOpt = tvnormOptions('AbsTol',1e-2);
tvOpt  = tvodeOptions;

%% Main For Loop
for i = 1:k
    % Time-varying representation
    Gss = Gall(:,:,i);
    G = evalt(tvss(Gss),[T0,Tf]);
    
    % Order
    Nx = order(Gss);
    
    % Power Iterations
    [pGain(i),~,info1]  = tvnormPower(G,NE,d1,pOpt);
    pTotalTime(i)       = info1.TotalTime;
    pTotalIter(i)       = info1.TotalIter;
    pTimeAllIter{i}     = info1.allIterComputeTime;
    pTemp               = pTimeAllIter{i};
    pSinglePowerIt(i)   = sum(pTemp(1:end-1))/length(pTemp(1:end-1));
    
    % RDE Bisection
    [tvn,~,info2]       = tvnormBisect(G,NE,tvnOpt);
    rUB(i)              = tvn(2);
    rLB(i)              = tvn(1);
    rTotalBisections(i) = info2.TotalBisections;
    rTotalTime(i)       = info2.TotalTime;
    
    % Use Combined Method
    [gc,~,info3]        = tvnormCombined(G,NE,d1,tvnOpt);
    cUB(i)              = gc(2);
    cLB(i)              = gc(1);
    cRDECount(i)        = info3.RDEcnt;
    cTotalTime(i)       = info3.TotalTime;
    cPowerIterations{i} = info3.Lower.PowerIterInfo;
    
    % Perform single CDRE integration
    gTry = rUB(i);
    [gSingle,dwc,info4] = tvnorm2(G,gTry,NE,tvOpt);
    
    % Make sure dwc is empty and gTry is indeed an upper bound
    if isempty(dwc)
        rSingleRDE(i) = info4.IntegrationTime;
    else
        warning('Single CDRE iteration did not converge.');
    end
    
    % Display iteration number
    disp(i);
end

% Save Data
filename = [mfilename sprintf('_nPlant%d',nPlant)];
save(filename,'pGain','pTotalTime','pTotalIter','pTimeAllIter','pSinglePowerIt',...
    'rUB','rLB','rTotalBisections','rTotalTime','rSingleRDE',...
    'cUB','cLB','cRDECount','cPowerIterations','cTotalTime');
return;

%% Plot Data
load('CompareTotalCompTime_nPlant5');
load('LTIDataSet_nPlant5','nPlant','NxAll');

% Average Per Model Order
k = 1;
for i = 1:nPlant:length(rTotalTime)
    rTemp = rTotalTime(i:i+nPlant-1);
    cTemp = cTotalTime(i:i+nPlant-1);
    pTemp = pTotalTime(i:i+nPlant-1);
    
    tB_AVG(k)  = sum(rTemp)/nPlant;
    tP_AVG(k)  = sum(pTemp)/nPlant;
    tC_AVG(k)  = sum(cTemp)/nPlant;
    
    tP_AVG_PerIter(k)  = sum(pSinglePowerIt(i:i+nPlant-1))/nPlant;
    tSingleCDRE_AVG(k) = sum(rSingleRDE(i:i+nPlant-1))/nPlant;
    k = k + 1;
end

% Plot Total Computational Time
figure(1);clf;
plot_fcn = @semilogy;
plot_fcn(NxAll,tP_AVG,'gd',NxAll,tB_AVG,'bo',NxAll,tC_AVG,'rs','LineWidth',2);

NxArray = 1:1:200;
f1 = fit(NxAll',log(tP_AVG'),'poly1'); hold on;
z1 = f1(NxArray');
p1 = plot_fcn(NxArray,exp(z1),'g');
p1.LineWidth = 2;

f2 = fit(NxAll',log(tB_AVG'),'poly1'); hold on;
z2 = f2(NxArray');
p2 = plot_fcn(NxArray,exp(z2),'b');
p2.LineWidth = 2;

f3 = fit(NxAll',log(tC_AVG'),'poly1'); hold on;
z3 = f3(NxArray');
p3 = plot_fcn(NxArray,exp(z3),'r');
p3.LineWidth = 2;

grid on;box on;
xlabel('System Order (n_x)');
ylabel('Average Computational Time (sec)');
%legend([p1,p2,p3],{'Power Iteration','RDE Bisection','Combined Algorithm'},'Location','southeast');
legend('Power Iteration','RDE Bisection','Combined Algorithm','Location','southeast');

% Plot ratio of RDE time/Power iter time
figure(2);clf;
tRatio = tSingleCDRE_AVG./tP_AVG_PerIter;
pr = plot(NxAll,tRatio,'b^','LineWidth',2);
fr = fit(NxAll',tRatio','poly1'); hold on;
pfr = plot(fr,'r');
pfr.LineWidth = 2;
grid on;box on;
legend([pr,pfr],{'Data','Linear fit'},'Location','southeast')
ylabel('T_{RDE}(n_x) / T_{PI}(n_x)');
xlabel('System Order (n_x)');