%% Load Data
load('FHvsIHSySingleStateEx.mat');
tvropt.MaxIter = 10;
GuncSCL = blkdiag(sqrt(beta),eye(NY-Nv))*Gunc*blkdiag(sqrt(beta),eye(NU-Nw));

%% Horizon (sec.)
TAll = [1,3,10,50,70,100,200,250,300];
NT = length(Tall);

% Memory Allocation
Krob  = cell(NT,1);
CLrob = cell(NT,1);
grob  = zeros(NT,1);

%% Finite Horizon L2toL2 Synthesis
for i = 1:NT
    Gtv = tvss(GuncSCL,[0 TAll(i)]);
    [Krob{i},CLrob{i},grob(i)] = tvrobsyn(Gtv,Delta,Ny,Nu,NE,tvropt);
end

%% Save
save(mfilename,'TAll','Krob','CLrob','grob');