%% TVNORM Example

% Plant P
load('RSSEx.mat','Pall');
P = Pall(:,:,1);

[A,B,C,D] = ssdata(P);
A = round(A*10)/10;
B = round(B*10)/10;
C = round(C*10)/10;
P = ss(A,B,C,0);
P = P(1,1);

tvnopt = tvnormOptions('OdeSolver','ode23s','Display','off');

%% Finite Horizon TVNORM
T0 = 0;
Ts = 0.1;
Tf = [0.5 1 2.5 5 8 12 20 30];
NE = 0;

tvnfh = zeros(size(Tf));
for i = 1:length(Tf)
    Gfh = evalt(tvss(P),T0:Ts:Tf(i));
    tvn = tvnorm(Gfh,NE,tvnopt);
    tvnfh(i) = tvn(2);
    fprintf(' [%s] Tf = %1.1f, Gain = %.4f\n\n',datetime,Tf(i),tvnfh(i))
end

%% Plot
figure(1);clf;hold on;box on;grid on;
plot(Tf,tvnfh,'-*b','LineWidth',2.5,'MarkerSize',8);
ylabel('Performance Metric','FontSize',14);
xlabel('Time (s)','FontSize',14);

if isequal(NE,0)
    tvnih = hinfnorm(P);
    fprintf(' [%s] Tf = Inf, Gain = %.4f\n',tvnih)
    plot(Tf,tvnih*ones(size(Tf)),'--r','LineWidth',2.5);
    legend('Finite Horizon','Infinite Horizon','Location','best','FontSize',14);
else
    legend('Finite Horizon','Location','best','FontSize',14);
end