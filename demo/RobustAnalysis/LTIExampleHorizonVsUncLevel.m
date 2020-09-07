%% LTIExampleHorizonVsUncLevel
% This example computes the worst-case gain for different horizons and
% different uncertainty level. At the end a 3D surface plot is obtained
% which provides interesting insights on how the worst-case gain varies as
% horizon and uncertainty level changes.

%% First order LTI system
% Create some lightly damped system
P = ss(0,1,1,0);
K = ss(1);
systemnames = 'P K';
inputvar = '[w; r]';
outputvar = '[r-K; P]';
input_to_P = '[w+r-K]';
input_to_K = '[P]';
Gs = sysic;
G = tvss(Gs);
Delta = ultidyn('Delta',[1 1]);
DeltaIQC = udyn('Delta',[1 1],'UserData',[0,-10,1]);

% Perform worst-case gain analysis for different horizons and different
% uncertainty level
T0 = 0;
Tall = [0.01 0.1 0.2 0.3 0.4 0.5 1   2   3   5   10];
Uall = [0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

NU = length(Uall);
NT = length(Tall);

% Options
tvwcopt = tvwcOptions;
wcopt = wcOptions;

% Allocate memory
g = cell(NU,NT);
gIH = cell(NU,1);
info = cell(NU,NT);

%% For loops
parfor i = 1:NU
    for j = 1:NT
        fprintf(' Running UL = %.1f, Tf = %.1f\n',Uall(i),Tall(j)); %#ok<PFBNS>
        Gj = evalt(G,[T0 Tall(j)]);
        
        % Uncertainty Level
        tvwcopt1 = tvwcopt;
        tvwcopt1.ULevel = Uall(i);
        
        [g{i,j},info{i,j}] = tvwcgain(Gj,DeltaIQC,0,tvwcopt1);
    end
    wcopt1 = wcopt;
    wcopt1.ULevel = Uall(i);
    gIH{i} = wcgain(lft(Delta,Gs),wcopt1);
end

%% Save Data
save(mfilename,'g','gIH','info','Tall','Uall');

%% Load Data & Plot
load('LTIExampleHorizonVsUncLevel.mat');
gmat = cell2mat(g);

% Presets
IDs = [3 6 9 10];
Marker = {'b*-','r^-','gs-','md-'};
legendArray = cell(1,length(IDs));

figure;hold on;
k = 1;
for i = IDs
    plot(Tall,gmat(i,:),Marker{k});
    legendArray{k} = ['\beta = ' sprintf('%.1f',Uall(i))];
    k = k + 1;
end
ylabel('Worst-case L_2 gain');
xlabel('Horizon (T)');
grid on; box on;hold off;
legend(legendArray,'location','southeast');

figure;hold on;k = 1;
for i = IDs
    plot(Uall,gmat(:,i),Marker{k});
    legendArray{k} = ['T = ' sprintf('%.1f',Tall(i))];
    k = k + 1;
end
ylabel('Worst-case L_2 gain');
xlabel('Unvcertainty Level (\beta)');
grid on; box on;hold off;
legend(legendArray,'location','northwest');

figure;
surface(Uall,Tall,gmat);
grid on;box on;colorbar;
xlabel('Unvcertainty Level (\beta)');
ylabel('Horizon (T)');
zlabel('Worst-case L_2 gain');