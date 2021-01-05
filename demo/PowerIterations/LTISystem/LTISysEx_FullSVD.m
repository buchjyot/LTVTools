%% LTISysEx_FullSVD
% This example computes all singular functions for LTI/LTV plant on finite
% horizons

% Seed
rng(0);

% Random System
Nx = 3;
Ne = 2;
Nd = 2;
Gss = rss(Nx,Ne,Nd);
Gss.D = 0;

% Horizon
T0 = 0;
Tf = 3;
Ts = 0.01;
t  = T0:Ts:Tf;
Nt = length(t);

% Time-varying system and adjoint
G = tvss(Gss,[T0,Tf]);
Ga = tvss(Gss',[T0,Tf]);
x0 = zeros(Nx,1);
lamTf = zeros(Nx,1);

% Options
Display = true;
StopTol = 5e-3;
MaxIter = 100;
tvopt = tvodeOptions('OdeSolver','ode45');

% Memory Allocation
d    = cell(MaxIter,1);
cost = zeros(MaxIter,1);

% Assume you want to compute first imax number of eigenvalues
imax  = 5;
psi   = cell(imax,1);
l2n   = zeros(imax,1);

% Induced L2 gain power iterations
for i = 1:imax
    
    % Display
    if Display
        fprintf('======================================\n');
        fprintf(' Approximating L2Gain: %d\n',i);
        fprintf('======================================\n');
    end
    
    % Choose Starting Vector
    d1 = tvmat(randn(Nd,1,Nt),t);
    for k = 1:i-1
        [psik,d1] = evalt(psi{k},d1,union(psi{k}.Time,d1.Time));
        alphak = psik'*d1;
        d1 = d1 - trapz(alphak.Time,alphak.Data)*psik; % Gram–Schmidt Process
    end
    d{1} = d1/tvnorm(d1);
    
    % Power iteration for loop for ith singular value
    for j = 1:MaxIter
        
        % State Equation
        e = tvlsim(G,d{j},[T0 Tf],x0,tvopt);
        
        % Costate Equation
        r = tvlsim(Ga,e,[Tf T0],lamTf,tvopt);
        
        % Compute Gain
        cost(j) = tvnorm(e); % Induced L2 gain, NOTE: d{j} is already unit norm
        
        % Display
        if Display
            fprintf(' Iter: %d, FwdPerf:%.3f\n',j,cost(j));
        end
        
        % Alignment Condition
        for k = 1:i-1
            [psik,r] = evalt(psi{k},r,union(psi{k}.Time,r.Time));
            alphak = psik'*r;
            r = r - trapz(alphak.Time,alphak.Data)*psik; % Gram–Schmidt Process
        end
        d{j+1} = r/tvnorm(r);
        
        % Stopping Condition (If input vector is stationary then stop)
        if (j > 1) && (j < MaxIter-1)
            if abs(cost(j)-cost(j-1)) <= StopTol
                [dj,djm1] = evalt(d{j},d{j-1},union(d{j}.Time,d{j-1}.Time));
                if tvnorm(dj-djm1) <= StopTol*2
                    psi{i} = d{j};
                    l2n(i) = cost(j);
                    break;
                end
            end
        end
    end
    
    % If the gain is approaching zero then stop, this will be the case only
    % if G is compact operator, i.e. zero feedthrough.
    if l2n(i)<=1e-6
        break;
    end
end

% Diagonal matrix of induced L2 gains, eigenvectors are psi
diag(l2n)