%% Power Iteration Example
% This example constructs a worst-case L2 disturbance using power iteration
% method

% Horizon
T0 = 0;
Tf = 3;

% Problem data
Nx = 2;
Ny = 1;
Nu = 1;
G = ss(tf(1,[1 0.6 1]));
G = ss(G.A,G.B,eye(2),0);

% Time-varying representation
Gt = evalt(tvss(G),[T0,Tf]);

% Options
NE = 2;
pOpt = poweritOptions('Display','on');
pSpec = poweritSignalSpec('InputL2Norm',1,'NE',NE);

%% Power Itertaions with zero Initial Condition
fprintf('### MATLAB Power Iterations with zero Initial Condition:\n');
[gLB1,dwc1,info1] = powerit(Gt,[T0,Tf],pSpec,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB1,info1.TotalIter,info1.TotalTime);

% Plot trajectories in phase plane
figure(1);clf;hold on;box on;grid on;

% trajac1
xwc1  = info1.Xwc;
xwc10 = tvsubs(xwc1,T0);
xwc1f = tvsubs(xwc1,Tf);
tvplot(xwc1(1),xwc1(2),'g',xwc1f(1),xwc1f(2),'go','LineWidth',2,'MarkerFaceColor','g');
plot(xwc10(1),xwc10(2),'ko','MarkerFaceColor','k');

% patch
[x1Ball1,x2Ball1] = getCircleData(xwc10,gLB1);
patch(x1Ball1,x2Ball1,'g');
alpha(0.1);

% Touch up the figure
axis equal;xlabel('x_1'); ylabel('x_2');
xlim([-1.5 1.5]);ylim([-1.5 1.5]);

%% Power Itertaions with non-zero Initial Condition
fprintf('### MATLAB Power Iterations with non-zero Initial Condition:\n');
x0 = [1.5; 1.5];
[gLB2,dwc2,info2] = powerit(Gt,[T0,Tf],x0,pSpec,pOpt);

% Display
fprintf(' Lower Bound: %.3f, Total Iterations: %d, Computational Time: %.3f seconds\n\n',gLB2,info2.TotalIter,info2.TotalTime);

% Plot trajectories in phase plane
figure(2);clf;hold on;box on;grid on;

% trajac2
xwc2  = info2.Xwc;
xwc20 = tvsubs(xwc2,T0);
xwc2f = tvsubs(xwc2,Tf);
tvplot(xwc2(1),xwc2(2),'r',xwc2f(1),xwc2f(2),'ro','LineWidth',2,'MarkerFaceColor','r');
plot(xwc20(1),xwc20(2),'ko','MarkerFaceColor','k');

% patch
[x1Ball2,x2Ball2] = getCircleData(xwc10,gLB2);
patch(x1Ball2,x2Ball2,'r');
alpha(0.1);

% Touch up the figure
axis equal;xlabel('x_1'); ylabel('x_2');
xlim([-2.5 2.5]);ylim([-2.5 2.5]);

%% Power iteration with non-zero Initial Condition and ||yF-y0|| terminal cost
keyboard; % NOTE: change the solver in powerit.m file, solver_case = 3
[gLB3,dwc3,info3] = powerit(Gt,[T0,Tf],x0,pSpec,pOpt);

% Plot trajectories in phase plane
figure(3);clf;hold on;box on;grid on;

% trajac1
xwc3  = info3.Xwc;
xwc30 = tvsubs(xwc3,T0);
xwc3f = tvsubs(xwc3,Tf);
tvplot(xwc3(1),xwc3(2),'c',xwc3f(1),xwc3f(2),'co','LineWidth',2,'MarkerFaceColor','c');
plot(xwc30(1),xwc30(2),'ko','MarkerFaceColor','k');

% patch
[x1Ball3,x2Ball3] = getCircleData(xwc30,gLB3);
patch(x1Ball3,x2Ball3,'c');
alpha(0.1);

% Touch up the figure
axis equal;xlabel('x_1'); ylabel('x_2');

%% Plot
% Disturbances
figure(4);clf;hold on;box on;
fh1 = tvplot(dwc1,'g'  ,'LineWidth',2.5);
fh2 = tvplot(dwc2,'r','LineWidth',2.5);
fh3 = tvplot(dwc3,'c','LineWidth',2.5);
xlabel('Time (sec)');grid on;
ylabel('Worst-Case Disturbance Signal (dwc)');
legend([fh1(1) fh2(1) fh3(1)],{'Opt1','Opt2','Opt3'},'Location','best');

% Plot all ellipses
% Plot trajectories in phase plane
figure(5);clf;hold on;box on;grid on;

% trajac1
tvplot(xwc1(1),xwc1(2),'g',xwc1f(1),xwc1f(2),'go','LineWidth',2,'MarkerFaceColor','g');
plot(xwc10(1),xwc10(2),'ko','MarkerFaceColor','k');

% patch
patch(x1Ball1,x2Ball1,'g');
alpha(0.1);

% trajac2
tvplot(xwc2(1),xwc2(2),'r',xwc2f(1),xwc2f(2),'ro','LineWidth',2,'MarkerFaceColor','r');
plot(xwc20(1),xwc20(2),'ko','MarkerFaceColor','k');

% patch
patch(x1Ball2,x2Ball2,'r');
alpha(0.1);

% trajac3
tvplot(xwc3(1),xwc3(2),'c',xwc3f(1),xwc3f(2),'co','LineWidth',2,'MarkerFaceColor','c');
plot(xwc30(1),xwc30(2),'ko','MarkerFaceColor','k');

% patch
patch(x1Ball3,x2Ball3,'c');
alpha(0.1);

% Touch up the figure
axis equal;xlabel('x_1'); ylabel('x_2');