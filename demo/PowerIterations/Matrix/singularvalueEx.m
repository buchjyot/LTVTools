%% Define Matrix M
rng(0);
matrix_case = 1;

switch matrix_case
    case 1
        % Take any random matrix
        Nr = 2;
        Nc = 3;
        M = randn(Nr,Nc);
        
    case 2
        % Identical singular values
        Nr = 4;
        Nc = 7;
        D = randn(3);
        [U,S,V] = svd(D);
        S(2,2) = S(1,1);
        M = U*S*V';
        
    case 3
        M = [1 0 0;0 1 0;0 0 1];
end

%% Power iteartions to compute maximum singular value
[sig1,u1s,v1s,info1] = svpowerit(M);

%% Singular values of M are the square root of the eigenvalues of M'*M or M*M'
% [lamMTM,v1,info2] = eigpowerit(M'*M);
% [lamMMT,u1,info3] = eigpowerit(M*M');
% disp(sqrt(lamMTM))
return;

%% Example for oscillatory (Eigen Value)
rng(0);
Niter = 25;
M = [1 1;0 -0.5]; % M = [1 0 0; 0 1 0; 0 0 0.98];
v = randn(2,1);
v = v/norm(v);
l = zeros(Niter,1);
for i = 1:Niter
    w = M*v;
    l(i) = norm(w);
    v = w/l(i);    
    if i > 1
        if abs(l(i) - l(i-1)) <= 1e-4
            break;
        end
    end
end
figure;
plot(1:i,l(1:i),'-ob','LineWidth',2)
grid on;box on;set(gca,'GridLineStyle','--');
xlabel('Iteration Count');
ylabel('Gain (s)');xlim([1 i]);
title('EigenValue Power Iteration')

%% Example for monotonic convergence (Singular Value)
rng(0);
Niter = 50;
M = randn(2,3); 
v = randn(3,1); 
v = v/norm(v);
s = zeros(Niter,1);
for i = 1:Niter
    w = M'*M*v;
    s(i) = norm(w);
    v = w/s(i);
    if i > 1
        if s(i) - s(i-1) <= 1e-4
            break;
        end
    end
end
figure;
plot(1:i,sqrt(s(1:i)),'-ob','LineWidth',2)
grid on;box on;set(gca,'GridLineStyle','--');
xlabel('Iteration Count');
ylabel('s^{1/2}');xlim([1 i]);
title('SVD Power Iteration Convergence')