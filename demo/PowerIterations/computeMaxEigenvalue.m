%% Power Iteration Example
% This example computes maximum eigenvalue using power iterations

% Dimention of matrix
m = 4;

% Square Matrix
M = rand(m);

% Maximum iterations
Niter = 100;

% Memory Allocation
z = zeros(m,Niter);
lam = zeros(Niter,1);

% Random vector
z(:,1) = rand(m,1);

% Power Iterations
for i = 1:Niter
    
    % State Equation
    z(:,i+1) = M*z(:,i);
    
    % Alignment condition
    lam(i) = norm(z(:,i+1));
    z(:,i+1) = z(:,i+1)/lam(i);
    
    % Stopping Condition
    if (i > 1) && (i < Niter-1)
        if abs(lam(i)-lam(i-1)) <= 1e-3
            break;
        end
    end
end

% Display
fprintf('MATLAB: %.3f, PowerIter: %.3f, TotalIter: %d\n',max(eig(M)),lam(i),i);
