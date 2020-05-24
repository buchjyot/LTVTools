%% Model
% Note: This model is configured to stop at T = 10 seconds.
model = 'NLSys';
load_system(model);

%% Power Iterations
rng(0);
pOpt = slpoweritOptions('InputL2Norm',3,'Display','on','MaxIter',50);
[glb,dwc,info] = slpowerit(model,pOpt);

%% Plot
figure;
tvplot(dwc,'LineWidth',2);