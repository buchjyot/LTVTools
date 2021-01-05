%% LTI Example
rng(1352);
Ts = 0.1;
plant_case = 2;
switch plant_case
    case 1
        Gss = ss(0.8,1,1,0,Ts);
    case 2
        Nx = 2; Ny = 2; Nu = 2;
        A = randn(Nx);A = A/max(svd(A));
        B = randn(Nx,Nu);
        C = randn(Ny,Nx);
        Gss = ss(A,B,C,0,Ts);
    case 3
        Gss = ss(1.1,1,1,0,Ts); % Unstable case
end
gInf = hinfnorm(Gss);

%% Horizon
T0 = 0;
Tf = [0.2, 0.5, 1, 3, 5, 10, 20];% 50, 100, 200];
NT = length(Tf);
dwc = cell(NT,1);
g = zeros(NT,1);
pOpt = poweritOptions('Display','on','StopTol',1e-4);
for i = 1:NT
    fprintf(' Tf = %.1f\n',Tf(i));
    Gtv = tvss(Gss,[T0 Tf(i)]);
    
    [g(i),dwc{i},info] = tvnormp(Gtv,0,[],[],pOpt);
end

%% Plot Disturbance

% Gain vs Horizon
figure,plot(Tf,g,'-bo',Tf,gInf*ones(NT,1),'--r','LineWidth',2);grid on;box on
xlabel('Horizon (sec)');
ylabel('Induced L_2 Gain');
legend('Finite Horizon','Infinite Horizon','Location','southeast');
xlim([Tf(1),Tf(end)]);

% Worst-case disturbance
distwc = dwc{end-1};
figure,tvstem(distwc,'LineWidth',2);grid on; box on;
ylabel('Worst-case disturbance');