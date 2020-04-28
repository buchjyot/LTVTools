%% Uncertain LTV Example:
% Compare various methods on finite and infinite horizon

%% Uncertain System
% Random Example
load LTIUncExData1;
G = MBad(:,:,7);
[Ag,Bg,Cg,Dg] = ssdata(G);

Ag = round(Ag*10)/10;
Bg = round(Bg*10)/10;
Cg = round(Cg*10)/10;
Dg = round(Dg*10)/10;
G = ss(Ag,Bg,Cg,Dg);

% Options
tvwcopt = tvwcOptions('Display','on');
Delta1 = udyn('Delta1',[1 1],'UserData',[0,0,0]);
Delta2 = udyn('Delta1',[1 1],'UserData',[1,0,0]);
NE = 0;

%% Robust L2 to L2 Gain on Finite Horizon
% Iterate between finite horizon LMI and RDE Solutions

% Perfom Analysis over a range of horizons
Tall = [1 2 5 10:10:50 100];
NT = numel(Tall);

% Analysis
Nall = zeros(NT,1);
RobTime = zeros(NT,1);
gfinal1 = zeros(NT,1);
gfinal2 = zeros(NT,1);
wcinfo1 = cell(NT,1);
wcinfo2 = cell(NT,1);

tstart = tic;
for k=1:NT
    % Finite Horizon
    T = Tall(k);
    fprintf('\n Analysis for T = %4.3f \t (%d of %d)\n',T,k,NT);
    
    % Robust L2 to L2 gain
    Gfh = evalt(tvss(G),linspace(0,T,tvwcopt.Nlmi));
    [gfinal1(k),wcinfo1{k}] = tvwcgain(Gfh,Delta1,NE,tvwcopt);
    [gfinal2(k),wcinfo2{k}] = tvwcgain(Gfh,Delta2,NE,tvwcopt);
end
tstop = toc(tstart);
save(mfilename,'gfinal1','gfinal2','wcinfo1','wcinfo2','Tall');

%% Plot Data
load('NLTVUncExample.mat','gfinal1','gfinal2','Tall')
plot(Tall,gfinal1,'--or',Tall,gfinal2,'-b^','LineWidth',2);
xlabel('Horizon (sec)','FontSize',14);title('Worst-case Induced L_2 Gain Upper Bound','FontSize',14);
ylabel('Induced L_2 Gain','FontSize',14);xlim([Tall(1) Tall(end-1)]);
legend('$\Delta \in \mathcal{S}_{TV}$','$\Delta \in \mathcal{S}_{ML}$','Location','southeast','interpreter','latex','FontSize',14);
grid on; box on;