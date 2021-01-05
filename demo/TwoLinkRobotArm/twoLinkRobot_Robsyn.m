%% twoLinkRobot_DKsyn
% Design a DK synthesis based controller for twoLinkRobotArm

%% Load twoLinkRobot_BuildLTVModel.mat file
load('twoLinkRobot_BuildLTVModel.mat');
load('twoLinkRobot_HinfDesign.mat');
load('twoLinkRobot_SpecifyOptions.mat');
[NY,NU] = size(Gunc);

%% Robust Synthesis with DelNorm = 0.8
% Scale Uncertainty Channel
Lscl = diag([sqrt(DelNorm) ones(1,NY-Nv)]);
Rscl = diag([sqrt(DelNorm) ones(1,NU-Nw)]);
Gunci = Lscl*Gunc*Rscl;

% Finite Horizon DK Synthesis Problem
[Krob1,CLrob1,grob1,robinfo1] = tvrobsyn(Gunci,Delta,Ny,Nu,NE,tvropt);

%% Robust Synthesis with DelNorm = 0.4
% Scale Uncertainty Channel
DelNorm = 0.4;
Lscl = diag([sqrt(DelNorm) ones(1,NY-Nv)]);
Rscl = diag([sqrt(DelNorm) ones(1,NU-Nw)]);
Gunci = Lscl*Gunc*Rscl;

% Finite Horizon DK Synthesis Problem
[Krob2,CLrob2,grob2,robinfo2] = tvrobsyn(Gunci,Delta,Ny,Nu,NE,tvropt);

%% Save Data
save(mfilename,'Gunc','Krob1','grob1','robinfo1','Krob2','grob2','robinfo2','-v7.3');