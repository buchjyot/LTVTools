%% twoLinkRobot_SpecifyOptions.m
% Specify options for ode tolerances, solver, display etc.
Display     = 'off';
OdeSolver   = 'ode23s';
RelTol      = 5e-3;
AbsTol      = 5e-4;
Bounds      = [0 100];

% General options
tvhopt = tvhinfsynOptions('Bounds',[0 5],'Display',Display,'OdeSolver',...
    OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvnopt = tvnormOptions('Bounds',Bounds,'Display',Display,'OdeSolver',....
    OdeSolver,'RelTol',RelTol,'AbsTol',AbsTol);
tvopt = tvodeOptions('OdeSolver',OdeSolver);
tvsopt = tvlsimOptions('OdeSolver',OdeSolver);
tvwcopt = tvwcOptions('RDEOptions',tvnopt,'Display','on','Nsp',10,'Nlmi',20,'StopTol',5e-3,'MaxIter',15);
tvropt = tvrobsynOptions('SynthesisOptions',tvhopt,'AnalysisOptions',tvwcopt,...
    'Display','on','MaxIter',7,'DebugMode',false);

%% Save
save(mfilename,'tvhopt','tvnopt','tvopt','tvwcopt','tvropt');