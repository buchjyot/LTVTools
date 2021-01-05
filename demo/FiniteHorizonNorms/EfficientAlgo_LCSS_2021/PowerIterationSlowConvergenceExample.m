%% Systems
% Horizon
T0 = 0;
Tf = 10;

% Seed
rng(0); %4545
MCase = 4;
switch MCase
    case 1
        M = eye(2);
        
    case 2
        [U,S,V] = svd(randn(2));
        S(2,2) = S(1,1);
        M = U*S*V';
        
    case 3
        [M,~,~] = svd(randn(2));
        
    case 4
        M = [1 0; 0 0.95];
end

% Systems
G1 = tf(1,[1 0.1 0.2]); hinfnorm(G1);
G2 = G1*M;

G1 = tvss(G1,[T0,Tf]);
G2 = tvss(G2,[T0,Tf]);

%% Power Iterations
pOpt = poweritOptions('Display','off','StoreAllIter',true,'StopTol',5e-3);
t = linspace(T0,Tf,100);
d1 = tvmat(ones(1,1,100),t);
[g1,dwc1,info1] = tvnormPower(G1,0,d1,[],pOpt);
[g2,dwc2,info2] = tvnormPower(G2,0,[d1;d1],[],pOpt);

%% RDE Approach
tOpt = tvnormOptions('Display','off','AbsTol',5e-3);
[n2,d2,rdeinfo1] = tvnormBisect(G2,tOpt);

%% Compare MIMO plant results to SISO plant results
[n1,d1,rdeinfo2] = tvnormBisect(G1,tOpt);
save(mfilename);
return;

%% Plot
figure(1);
niter1 = length(info1.allPerf);
niter2 = length(info2.allPerf);
plot(1:niter1,info1.allPerf,'-ob',1:niter2,info2.allPerf,'--sm','LineWidth',2);
grid on;box on;
xlabel('Iteration Count (i)');
ylabel('Induced L_2 Gain');
% ylabel('Induced L_2 Gain (\gamma_f^{(i)})');
legend('G_1','G_2','Location','southeast');
xlim([1 niter2]);
return;

%% Extra
% linear sim gain
[lsimg1,e1] = tvlsimgain(G1,d1);

% Take random vector
v = randn(2,1);

% Linear sim (lsimg2 must be close to n2)
[lsimg2,e2] = tvlsimgain(G2,d1*v);

%% Plot d1 and d2
figure(1);
tvplot(d1,'b','LineWidth',2);
ylabel('d1');grid on, box on;

figure(2);
tvplot(d2,'LineWidth',2);
ylabel('d2');grid on, box on;
legend('d2_1','d2_2');

%% Guess the output
e2guess = e1*M*v;

figure(3);
tvplot(e2guess(1),'b',e2guess(2),'r',e2(1),'--g',e2(2),'--k','LineWidth',2);
legend('e2guess_1','e2guess_2','e2_1','e2_2');title('Outputs')
grid on, box on;