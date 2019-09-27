%% TVNORM Example

% Load the plant
load('RSSEx.mat','P','T0','Ts');
tvnopt = tvnormOptions('OdeSolver','ode23s');

%% Finite Horizon TVNORM
Tf = [0.5 1 2.5 5 8 12 20 30];
tvnfh = zeros(size(Tf));
parfor i = 1:length(Tf)
    Gfh = evalt(tvss(P),T0:Ts:Tf(i));
    tvn = tvnorm(Gfh,tvnopt);
    tvnfh(i) = tvn(2);
    fprintf('Tf=%1.1f, Gain=%.4f\n',Tf(i),tvnfh(i))
end

%% Inf Horizon Results
tvnih = hinfnorm(P);
fprintf('Tf=Inf, Gain=%.4f\n',tvnih)
figure;clf;hold on;box on;grid on;
plot(Tf,tvnfh,'-*b','LineWidth',2.5,'MarkerSize',8);
plot(Tf,tvnih*ones(size(Tf)),'--r','LineWidth',2.5);
ylabel('Performance Metric','FontSize',14);
xlabel('Time (s)','FontSize',14);
legend('Finite Horizon','Infinite Horizon','Location','best','FontSize',14);