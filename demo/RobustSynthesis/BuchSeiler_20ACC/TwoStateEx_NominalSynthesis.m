%% Load Data
load('TwoStateEx.mat');

%% Nominal Synthesis
fprintf('====================\n');
fprintf(' Nominal Synthesis\n');
fprintf('====================\n');

% Design
[Knom,CLnom,gnom,nominfo] = tvhinfsyn(Gnom,Ny,Nu,NE,tvhopt);
fprintf(' Synthesis (Gamma):= %.4f\n',gnom);

% Verify Synthesis Results
[gn,dWcn] = tvnorm(CLnom,NE,tvnopt);
fprintf(' Analysis (Gamma):= %.4f\n',gn(2));

% Report computation time and coupled RDE bisections
fprintf(' CompTime:= %.4f, RDECnt:= %d\n',nominfo.TotalTime,nominfo.RDEcnt);

% Nominal Controller (With Design Weights)
CLn = lft(evalt(Gunc,Knom.Time),Knom);

%% Save Data
save(mfilename,'Knom','CLn','gnom','gn','dWcn','nominfo');