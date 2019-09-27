%% TVH2NORM Example
%
% This Example tests tvh2norm function which calculates finite horizon 2
% norm of a LTV system

% Load System Data
load('RSSEx.mat');

%% Finite Horizon 2-Norm
Tf = [0.5 1 3 5 8 15 20 30 50 70 100 200];
h2nfh = zeros(length(Tf),1);
for i = 1:length(Tf)
    Pi = evalt(tvss(P),T0:Ts:Tf(i));
    h2nfh(i) = tvh2norm(Pi,tvopt);
    fprintf(' Tf = %.1f, Closed Loop H2 Norm = %.4f\n',Tf(i),h2nfh(i));
end

%% Infinite Horizon 2-Norm
h2nih = norm(P,2);
fprintf(' Tf = Inf, Closed Loop H2 Norm = %.4f\n',h2nih);

%% Plot
if length(Tf) > 1
    figure;clf;
    plot(Tf,h2nih*ones(size(Tf)),'r--','LineWidth',3);hold on;box on;grid on;
    plot(Tf,h2nfh,'b*-','LineWidth',3,'MarkerSize',11);
    xlabel('Time Horizon (s)','FontSize',12);
    ylabel('H_2 Norm','FontSize',12);
    title('H_2 Norm on Finite Horizon','FontSize',12);
    legend('Infinite Horizon','Finite Horizon','Location','Best');
    XLIM = [0 Tf(end-1)];
    xlim(XLIM);
end