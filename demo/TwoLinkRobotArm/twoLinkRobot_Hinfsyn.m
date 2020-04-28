%% twoLinkRobot_Hinfsyn
% Design Hinfinity Output Feedback Controller  & perform analysis

%% Load twoLinkRobot_BuildLTVModel.mat file
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_SpecifyOptions.mat');
tvnopt.RelTol = 5e-3;

%% Design Output Feedback Hinf Controller

% Synthesis Step
[Knom,CLnom,gnom,nominfo] = tvhinfsyn(Gnom,Ny,Nu,NE,tvhopt);
fprintf(' Synthesis TVHINFSYN Bound:');
disp(gnom);

% Verify bisection results
[g,dwc] = tvnorm(CLnom,NE,tvnopt);
fprintf(' TVNORM Bisection Bound:');
disp(g(2));

%% Save Data
save(mfilename,'Knom','Gnom','gnom','nominfo');