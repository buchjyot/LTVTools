function [Y,X,U] = tvstep(G,Tfinal,Opt)
% TVSTEP simulate the step response of LTV System G
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
Tf = GTime(end);
switch nin
    case 1
        Tfinal = Tf;
        Opt = tvstepOptions;
    case 2
        Opt = tvstepOptions;
end

% Assume Zero Initial Conditions
Nx = order(G);
x0 = zeros(Nx,1);
[~,Nu] = size(G);

% Construct input U
U = LOCALGetU(GTime,Tfinal,Nu,Opt);

% Simulate Response
tvOpt = tvodeOptions('OdeSolver',Opt.OdeSolver,'OdeOptions',Opt.OdeOptions);
[Y,X] = tvlsim(G,U,x0,tvOpt);

function U = LOCALGetU(GTime,Tfinal,Nu,Opt)
% Tfinal must be between T0 and Tf and must be greater than Tstart
Tstart = Opt.StepTime;
T0 = GTime(1);
Tf = GTime(end);
if Tfinal > Tf || Tfinal < T0 || Tfinal < Tstart
    error('Invalid Tfinal');
end

UST = cat(1,Tstart,GTime(GTime>Tstart & GTime<Tfinal),Tfinal);
StepAmp = Opt.StepAmplitude;
U = tvmat(StepAmp*ones(Nu,length(UST)),UST);

warning('off','ltvtools:evalt:extrapolate');
% Sharpen the edge
UT = cat(1,U.Time(1)-0.003,U.Time(1)-0.002,U.Time(1)-0.001,...
    U.Time,U.Time(end)+0.001,U.Time(end)+0.002,U.Time(end)+0.003);
% Transform to GTime Grid
U = evalt(U,union(GTime,UT));
Uadd = tvmat(Opt.InputOffset*ones(Nu,length(U.Time)),U.Time);
U = U + Uadd;
warning('on','ltvtools:evalt:extrapolate');