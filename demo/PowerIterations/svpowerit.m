function [sig1,u1,v1,info] = svpowerit(M)
%% SVDPOWERIT
% This example computes maximum singularvalue and associated singular
% vectors using power iterations

% Size of M
[m,n] = size(M);

% Display
Display = true;

% Maximum iterations
MaxIter = 500;

% Stopping Tolerance
StopTol = 1e-2;
StopTolSatisfied = false;

%% Memory Allocation
v       = zeros(n,MaxIter);  % First Right Singular Vector
u       = zeros(m,MaxIter);  % First Left Singular Vector
sigma   = zeros(MaxIter,1);

% Random Vector (Normalized)
v(:,1) = rand(n,1);
v(:,1) = v(:,1)/norm(v(:,1));

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    u(:,i) = M*v(:,i);
    
    % Compute Gain
    sigma(i) = norm(u(:,i));
    
    % Alignment Condition
    u(:,i) = u(:,i)/ sigma(i);
    
    % Display
    if Display
        fprintf(' Iter: %d, Sigma1: %.3f\n',i,sigma(i));
    end
    
    % Stopping Condition (If input vector is stationary then stop)
    if (i > 1) && (i < MaxIter-1)
        if abs(sigma(i)-sigma(i-1)) <= StopTol
            if norm(v(:,i)-v(:,i-1))<= StopTol
                StopTolSatisfied = true;
                break;
            end
        end
    end
    
    % Costate Equation
    v(:,i+1) = M'*u(:,i);
    
    % Alignment Condition
    v(:,i+1) = v(:,i+1)/ norm(v(:,i+1));
    
end
tTotal = toc(t1);

%% Outputs
sig1 = sigma(i);
u1 = u(:,i);
v1 = v(:,i);
if Display
    if isequal(i,MaxIter) && ~StopTolSatisfied
        fprintf(' Maximum number of iterations reached.\n');
    end
    if StopTolSatisfied
        fprintf(' Stopping tolerance satisfied, terminating iterations.\n');
    end
end

info = [];
info.TotalIter = i;
info.TotalTime = tTotal;
info.allSigma = sigma(1:i);
end