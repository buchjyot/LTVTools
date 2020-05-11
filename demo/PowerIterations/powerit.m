function [glb,dwc,info] = powerit(G,varargin)
%% POWERIT
% Function that performs power iterations to compute worst-case gain lower
% bound for LTV system on a finite horizon
%
% Input:
% G     = must be TVSS on finite horizon
% Opt   = options for power iterations
%
% Output:
% glb   = lower bound on finite horizon induced L2 norm
% dwc   = final worst-case disturbance
% info  = additional info

%% Input Processing
narginchk(1,2);
nin = nargin;

% Read varargin
switch nin
    case 1
        Opt = poweritOptions;
    case 2
        if isa(varargin{1},'poweritOptions')
            Opt = varargin{1};
        else
            error('Options must be specified using "poweritOptions".');
        end
end

% Read Options
MaxIter         = Opt.MaxIter;
StopTol         = Opt.StopTol;
OdeSolver       = Opt.OdeSolver;
OdeOptions      = Opt.OdeOptions;
Display         = isequal(Opt.Display,'on');
StepSize        = Opt.StepSize;
tvsopt1         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','ForwardInTime','StepSize',StepSize);
tvsopt2         = tvlsimOptions('OdeSolver',OdeSolver,'OdeOptions',OdeOptions,'Type','BackwardInTime','StepSize',StepSize);
FixedStepSize   = isa(StepSize,'double');

% Use Plant Time Grid
Tgrid = G.Time;
Nt = length(Tgrid);

% System dimentions
[~,Nu] = size(G);

% Adjoint System
Ga = G';

%% Memory Allocation
U   = cell(MaxIter,1);
Y   = cell(MaxIter,1);
l2n = zeros(MaxIter,1);

% Select a Random Input Signal for first iteration
U{1} = tvmat(rand(Nu,1,Nt),Tgrid);
U{1} = U{1}/tvnorm(U{1});

%% Power Itertaions
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    Y{i} = tvlsim(G,U{i},tvsopt1);
    
    % Compute L2 Norm of the output
    l2n(i) = tvnorm(Y{i});
    
    % Display
    if Display
        fprintf('Iter: %d, l2norm: %.3f\n',i,l2n(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter-1)
        % Performance that we care is stationary
        if abs(l2n(i)-l2n(i-1)) <= StopTol*max(l2n([i i-1]))
            % Performance may reach stationary but input signal may not
            if FixedStepSize
                if tvnorm(U{i}-U{i-1}) <= StopTol
                    break;
                end
            else
                % For variable stepsize, we will have two different time
                % grids, thus for subtraction we need consistent time grid.
                tUnion = union(U{i}.Time,U{i-1}.Time);
                [Ui,Uim1] = evalt(U{i},U{i-1},tUnion);
                if tvnorm(Ui-Uim1) <= StopTol
                    break;
                end
            end
        end
    end
    
    % Alignment Condition
    Y{i} = Y{i}/l2n(i);
    
    % Costate Equation
    U{i+1} = tvlsim(Ga,Y{i},tvsopt2);
    
    % Alignment Condition
    U{i+1} = U{i+1}/tvnorm(U{i+1});
end
tcomp = toc(t1);

%% Final Output
glb = l2n(i);
dwc = U{i};

info = [];
info.TotalIter  = i;
info.TotalTime  = tcomp;
end