%% FullSVDEx
% This file computes full SVD of a random matrix

% Set seed
rng(0);

% Random matrix
nr = 3;
nc = 3;
M = randn(nr,nc);
mdim = min(nr,nc);

% Options
Display = true;
StopTol = 1e-4;
MaxIter = 100;

% Memroy allocation
v       = zeros(nc,MaxIter);  % First Right Singular Vector
sigma   = zeros(MaxIter,1);

% We will fill the following variables
V = zeros(nc);
S = zeros(nr,nc);

%% Singularvalue power iterations
for i = 1:mdim
    
    % Display
    if Display
        fprintf('======================================\n');
        fprintf(' Approximating Singular Value: %d\n',i);
        fprintf('======================================\n');
    end
    
    % Choose Starting Vector
    v1 = rand(nc,1);
    for k = 1:i-1
        v1 = v1 - (V(:,k)'*v1)*V(:,k); % Gram–Schmidt Process
    end
    v(:,1) = v1/norm(v1);
    
    % Power iteration for loop for ith singular value
    for j = 1:MaxIter
        
        % State & Costate Equations Together
        r = M'*M*v(:,j);
        
        % Compute Gain
        sigma(j) = sqrt(norm(r));
        
        % Display
        if Display
            fprintf(' Iter: %d, Sigma%d: %.3f\n',j,i,sigma(j));
        end
        
        % Alignment Condition
        for k = 1:i-1
            r = r - (V(:,k)'*r)*V(:,k); % Gram–Schmidt Process
        end
        v(:,j+1) = r/norm(r);
        
        % Stopping Condition (If input vector is stationary then stop)
        if (j > 1) && (j < MaxIter-1)
            if abs(sigma(j)-sigma(j-1)) <= StopTol
                if norm(v(:,j)-v(:,j-1))<= StopTol
                    V(:,i) = v(:,j);
                    S(i,i) = sigma(j);
                    break;
                end
            end
        end
    end
    
end

% Approximate Left Singular Matrix
U = M*V/S
S
V

%% Compare results with MATLAB's SVD
[mU,mS,mV] = svd(M)