function [lam1,v1,info] = eigpowerit(M)
%% EIGPOWERIT
% This function computes maximum eigenvalue and associated eigenvector
% using power iterations
%
% Input
% Square matrix M of size (m x m)
%
% Output
% lam1 : Maximum eigenvalue
% v1   : Associated eigenvector
% info : additional info

% Maximum iterations
MaxIter = 500;

% Display
Display = true;

% Stopping Tolerance
StopTol = 1e-2;
StopTolSatisfied = false;

%% Memory Allocation
m = size(M,1);
n = size(M,1);
if ~isequal(m,n)
    error('Input matrix must be square.');
end
z = zeros(m,MaxIter);
lam = zeros(MaxIter,1);

% Random vector
z(:,1) = rand(m,1);
z(:,1) = z(:,1)/norm(z(:,1));

%% Power Iterations
t1 = tic;
for i = 1:MaxIter
    
    % State Equation
    z(:,i+1) = M*z(:,i);
    
    % Alignment condition
    lam(i) = norm(z(:,i+1));
    z(:,i+1) = z(:,i+1)/lam(i);
    
    % Display
    if Display
        fprintf(' Iter: %d, Lam1: %.3f\n',i,lam(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter-1)
        if abs(lam(i)-lam(i-1)) <= StopTol
            if norm(z(:,i)-z(:,i-1)) <= StopTol
                StopTolSatisfied = true;
                break;
            end
        end
    end
end
tTotal = toc(t1);

%% Outputs
lam1 = lam(i);
v1 = z(:,i);
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
info.allLam = lam(1:i);
end