%% Random Constant Matrix
rng(52);
nd = 2;
M0 = randn(nd);
M0 = M0/norm(M0);

% Horizon
T0 = 0;
Tf = 5;
t  = T0:0.001:Tf;

% Random matrix
S = tvmat(sin(t),t);
C = tvmat(cos(t),t);
M = M0 + M0*blkdiag(S,C)/2; % M = tvmat(M0,t);
Nt = length(t);

% Random vector
v = tvmat(ones(nd,1,Nt),t);
v = v/tvnorm(v);

% Options and Memory Allocation
MaxIter = 100;
StopTol = 0.001;
fwcost = zeros(MaxIter,1);
bwcost = zeros(MaxIter,1);

% Power iteration
for i = 1:MaxIter
    % Forward
    vNext = M*v;
    fwcost(i) = tvnorm(vNext);
    
    % Backward
    uNext = M'*vNext/fwcost(i);
    bwcost(i) = tvnorm(uNext);
    
    % Update
    v = uNext/bwcost(i);
    
    figure(1);
    tvplot(v,'b');
    pause(2);
    
    % Stop Condition
    if i > 1
        if fwcost(i)-fwcost(i-1) < StopTol            
            break;
        end
    end
end

%% Perform SVD
[U,S,V] = svd(M);
[Sm,id] = max(S);
tSm = tvmax(Sm);
disp([fwcost(i) tSm])

figure;
plot(1:i, fwcost(1:i),'-ob',1:i,bwcost(1:i),'-rs');
xlabel('Iteration Count');
ylabel('Cost');
xlim([1 i]);