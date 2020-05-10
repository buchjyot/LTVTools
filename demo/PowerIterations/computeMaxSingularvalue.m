%% Power Iteration Example
% This example computes maximum singularvalue using power iterations

% Dimention of matrix
m = 3;
n = 5;

% Square Matrix
M = rand(m,n);
[U,S,V] = svd(M);

% Maximum iterations
Niter = 100;

% Memory Allocation
v       = zeros(n,Niter);  % First Right Singular Vector
u       = zeros(m,Niter);  % First Left Singular Vector
sigma   = zeros(Niter,1);

% Stopping Tolerance
StopTol = 1e-3;

% Random Vector (Normalized)
v(:,1) = rand(n,1);
v(:,1) = v(:,1)/norm(v(:,1));

% Power Iterations
for i = 1:Niter
    
    % State Equation
    u(:,i) = M*v(:,i);
    
    % Compute Gain
    sigma(i) = norm(u(:,i));
    
    % Alignment Condition
    u(:,i) = u(:,i)/ sigma(i);
    
    % Stopping Condition (If input vector is stationary then stop)
    if (i > 1) && (i < Niter-1)
        if norm(v(:,i)-v(:,i-1))<= StopTol
            break;
        end
    end
    
    % Costate Equation
    v(:,i+1) = M'*u(:,i);
    
    % Alignment Condition
    v(:,i+1) = v(:,i+1)/ norm(v(:,i+1));
    
end

% Display
fprintf('MATLAB: %.3f, PowerIter: %.3f, TotalIter: %d\n',max(max(S)),sigma(i),i);

% Singularvectors are close enough
fprintf('MATLAB v1 = \n');disp(V(:,1));
fprintf('PowerIter v1 = \n');disp(v(:,i));

fprintf('MATLAB u1 = \n');disp(U(:,1));
fprintf('PowerIter u1 = \n');disp(u(:,i));