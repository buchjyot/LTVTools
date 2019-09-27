%% TVGRAM Example
%
% This file tests Time-varying gramians for a random system

% Load System Data
load('RSSEx.mat','P','T0','Ts','Tf');

% GRAM
Wcih = gram(P,'c');
Woih = gram(P,'O');

% TVGRAM
tvgopt = tvodeOptions('OdeSolver','ode23s');
Pfh = evalt(tvss(P),T0:Ts:4);
Wcfh = tvgram(Pfh,'c',tvgopt);
Wofh = tvgram(Pfh,'O',tvgopt);

% Plotting
figure(1);clf;box on;grid on;hold on;
p1=tvplot(Wcfh,'b','LineWidth',2.5);
p2=plot(Wcfh.Time,repmat(Wcih(:),1,length(Wcfh.Time)),'--r','LineWidth',2.5);
title('Controllability Gramian');
legend([p1(1);p2(1)],'Finite Horizon','Infinite Horizon','Location','best');
xlabel('Time(s)');ylabel('W_c');

figure(2);clf;box on;grid on;hold on;
p1=tvplot(Wofh,'b','LineWidth',2.5);
p2=plot(Wofh.Time,repmat(Woih(:),1,length(Wofh.Time)),'--r','LineWidth',2.5);
title('Observability Gramian');
legend([p1(1);p2(1)],'Finite Horizon','Infinite Horizon','Location','best')
xlabel('Time(s)');ylabel('W_o');