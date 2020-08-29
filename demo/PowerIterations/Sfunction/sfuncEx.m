%% Horizon
T0 = 0;
Tf = 3;

%% Options
pOpt = poweritOptions('Display','on');
pSpec = poweritSignalSpec('NE',1);

%% Power Iterations
rng(0);
x0 = 0*randn(2,1);
[glb,dwc,info] = powerit(@LTISysSfunc,[T0,Tf],x0,pSpec,pOpt);