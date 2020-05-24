function [glb,dwc,info] = ltipowerit(G,Tgrid)
%% LTIPowerIT
% Power iteration algorithm for LTI systems (ss objects)
%
% Input:
% G : LTI System (ss object)
% Tgrid: Time grid for integration
%
% Output:
% glb: lower bound on induced l2 gain
% dwc: worst-case disturbance
% info: additional info
%
% This function does not use any of the LTVTools functions.

%% Init
MaxIter = 500;
StopTol = 1e-2;
Display = true;
OdeSolver = @ode45;
OdeOpt    = odeset('RelTol',1e-3,'AbsTol',1e-6);
StopTolSatisfied= false;

%% Dynamics
[A,B,C,D] = ssdata(G);
[Aa,Ba,Ca,Da] = ssdata(G');

% Sizes
Nt = length(Tgrid);
Nu = size(G,2);
Nx = order(G);

%% Perfromance
% Induced L2 gain
fhL2Norm = @(time,data) sqrt(trapz(time,LOCALNorm(data,length(time)).^2));

%% Memory Allocation
U   = cell(MaxIter,1);
Y   = cell(MaxIter,1);
l2n = zeros(MaxIter,1);

%% Initial Input
x0 = zeros(Nx,1);
lam0 = zeros(Nx,1);
U1 = randn(Nt,Nu);
U{1} = U1/fhL2Norm(Tgrid,U1);

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    % State Equation
    Uts = timeseries(U{i},Tgrid);
    [Tgrid,x] = OdeSolver(@(t,x) LOCALlsim(t,x,Uts,A,B),Tgrid,x0,OdeOpt);
    Ui = get(resample(Uts,Tgrid),'Data');
    Y{i} = (C*x' + D*Ui')';
    Yts = timeseries(Y{i},Tgrid);
    
    % Normalize
    l2n(i) = fhL2Norm(Tgrid,Y{i});
    Y{i} = Y{i}/l2n(i);
    
    % Display
    if Display
        fprintf(' Iter: %d, L2Gain: %.3f\n',i,l2n(i));
    end
    
    % Stop Condition
    if i > 1 && (i < MaxIter-1)
        % Check if performance that we care is stationary
        if abs(l2n(i)-l2n(i-1)) <= StopTol*max(l2n([i i-1]))
            if fhL2Norm(Tgrid,abs(U{i}-U{i-1})) <= StopTol
                StopTolSatisfied = true;
                break;
            end
        end
    end
    
    % Costate Equation
    AdjTime = flip(Tgrid);
    [Tgrid,lam] = OdeSolver(@(t,lam) LOCALlsim(t,lam,Yts,Aa,Ba),AdjTime,lam0,OdeOpt);
    Tgrid = flip(Tgrid);
    Yi = get(resample(Yts,Tgrid),'Data');
    lam = flipud(lam);
    U{i+1} = (Ca*lam' + Da*Yi')';
    
    % Normalize
    U{i+1} = U{i+1}/fhL2Norm(Tgrid,U{i+1});
end
tTotal = toc(t1);

%% Process Outputs
glb = l2n(i);
dwc = timeseries(U{i},Tgrid,'Name','Worst-Case Disturbance (dwc)');
if Display
    if isequal(i,MaxIter) && ~StopTolSatisfied
        fprintf(' Maximum number of iterations reached.\n');
    end
    if StopTolSatisfied
        fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
    end
end

info = [];
info.TotalTime = tTotal;
info.TotalIter = i;
info.allPerf = l2n(1:i);
end

function xdot = LOCALlsim(t,x,u,A,B)
uData = get(resample(u,t),'Data');
xdot = A*x + B*uData';
end

function out = LOCALNorm(V,Nt)
out = zeros(Nt,1);
for i = 1:Nt
    out(i) = norm(V(i,:),2);
end
end