%% Load Data
load('TwoStateEx.mat')

%% Robust Synthesis
fprintf('====================\n');
fprintf(' Robust Synthesis\n');
fprintf('====================\n');

% Weight Uncertainty channels using beta
[NY,NU] = size(Gunc);
GuncSCL = blkdiag(sqrt(beta),eye(NY-Nv))*Gunc*blkdiag(sqrt(beta),eye(NU-Nw));
GuncSCL.UserData = Gunc.UserData;

% Design
[Krob,CLrob,grob,robinfo] = tvrobsyn(GuncSCL,Ny,Nu,NE,tvdkopt);
fprintf('Synthesis (Gamma):= %.4f\n',grob);

% Robust Controller (With Design Weights)
CLr = lft(evalt(Gunc,Krob.Time),Krob);
CLr.UserData = Gunc.UserData;

%% Save Data
save(mfilename,'Krob','CLr','grob','robinfo');