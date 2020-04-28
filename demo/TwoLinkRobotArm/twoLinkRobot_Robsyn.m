%% twoLinkRobot_DKsyn
% Design a DK synthesis based controller for twoLinkRobotArm

%% Load twoLinkRobot_BuildLTVModel.mat file
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_SpecifyOptions.mat');

%% Robust Synthesis
% Scale Uncertainty Channel
[NY,NU] = size(Gunc);
Lscl = diag([sqrt(DelNorm) ones(1,NY-Nv)]);
Rscl = diag([sqrt(DelNorm) ones(1,NU-Nw)]);
Gunci = Lscl*Gunc*Rscl;

% Finite Horizon DK Synthesis Problem
[Krob,CLrob,grob,robinfo] = tvrobsyn(Gunci,Delta,Ny,Nu,NE,tvropt);

%% Save Data
save(mfilename,'Gunc','Krob','grob','robinfo','-v7.3');