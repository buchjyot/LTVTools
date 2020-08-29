function [Y,X] = tvstep(G,Tf,Opt)
%% TVSTEP simulate the step response of LTV System G
%
% Example:
%
% >> Gss = ss(-1,1,1,0);
% >> G = evalt(tvss(Gss),0:0.1:10);
% >> tsopt = tvstepOptions('InputOffset',1,'StepAmplitude',2,'StepTime',1);
% >> [y,x] = tvstep(G,5,tsopt)

% Input Processing
narginchk(1,3);
nin = nargin;
GTime = G.Time;
GT0 = GTime(1);
GTf = GTime(end);
switch nin
    case 1
        Tf = GTf;
        Opt = tvstepOptions;
    case 2
        Opt = tvstepOptions;
end
T0 = GT0;

% Assume Zero Initial Conditions
Nx = order(G);
x0 = zeros(Nx,1);
[Ny,Nu] = size(G);

% Memory Allocation
Y = cell(Ny,Nu);
X = cell(Ny,Nu);

% Construct input U
U = LOCALGetU(GTime,T0,Tf,Opt);

% Simulate Response
tvOpt = tvodeOptions('OdeSolver',Opt.OdeSolver,'OdeOptions',Opt.OdeOptions);

for i = 1:Ny
    for j = 1:Nu
        Gi = tvss(G.Data(i,j,:),GTime,G.InterpolationMethod);
        [Y{i,j},X{i,j}] = tvlsim(Gi,U,[T0,Tf],x0,tvOpt);        
    end
end
end

function U = LOCALGetU(GTime,T0,Tf,Opt)
% Tfinal must be between T0 and Tf and must be greater than Tstart
Tstart = Opt.StepTime;
if Tf > Tf || Tf < T0 || Tf < Tstart
    error('Invalid Tfinal');
end

UST = cat(1,Tstart,GTime(GTime>Tstart & GTime<Tf),Tf);
StepAmp = Opt.StepAmplitude;
U = tvmat(StepAmp*ones(1,length(UST)),UST);

warning('off','ltvtools:evalt:extrapolate');
% Sharpen the edge
UT = cat(1,U.Time(1)-0.003,U.Time(1)-0.002,U.Time(1)-0.001,...
    U.Time,U.Time(end)+0.001,U.Time(end)+0.002,U.Time(end)+0.003);
% Transform to GTime Grid
U = evalt(U,union(GTime,UT));
Uadd = tvmat(Opt.InputOffset*ones(1,length(U.Time)),U.Time);
U = U + Uadd;
warning('on','ltvtools:evalt:extrapolate');
end