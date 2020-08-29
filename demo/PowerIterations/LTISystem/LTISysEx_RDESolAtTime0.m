%% Example
% Horizon
T0 = 0;
Tf = 3;

% Problem data
G = ss(tf(1,[1 0.6 1]));
G = ss(G.A,G.B,eye(2),0);

% Sizes
Nx = 2;
Nu = 1;
NE = 2;

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

%% Compute TVNORM
tOpt = tvnormOptions('Display','on');
[n1,d1,info] = tvnorm(Gt,NE,tOpt);

%% Compte RDE Solutions

gGrid = [0.849 0.85 0.855 0.86 0.87 0.9 1];
Ng = length(gGrid);
P = cell(Ng,1);
P0 = cell(Ng,1);
figure(1);clf;grid on;box on;hold on;
for i = 1:Ng
   [~,~,info1] = tvnorm2(Gt,gGrid(i),NE,tOpt); 
   P{i} = info1.P;
   P0{i} = tvsubs(P{i},T0);
   
   p1 = plot_gaussian_ellipsoid([0,0],P0{i});
   p1.LineWidth = 2;
end
hold off;
xlabel('x_1');ylabel('x_2');
title(' RDE solution P(0) Ellipsoids')