%% TVMAT Interpolation Example
% This example shows you the impact of interpolation method for TVMAT using
% a simple example.

% Horizon
T0 = pi/4;
Tf = 7*pi/4;
f = @(t) sin(t);

%% True Data

% Pick a Very Dense Grid
TrueT = T0:0.00001:Tf;
TrueD = f(TrueT);
Tm = tvmat(TrueD,TrueT);

%% Sampled Data

% Create TVMAT using Corse Grid
ATime = T0:pi/4:Tf;
AData = f(ATime);
A = tvmat(AData,ATime);

% Pick a Dense Grid
TGrid = 0:pi/1600:2*pi;

%% Interpolation

% Set the warning to off
warning('off','ltvtools:evalt:extrapolate');

% Linear 
A.InterpolationMethod = 'Linear';
ALinIM = evalt(A,TGrid);

% Spline 
A.InterpolationMethod = 'Spline';
AsplIM = evalt(A,TGrid);

% Nearest
A.InterpolationMethod = 'Nearest';
ANeaIM = evalt(A,TGrid);

% Flat
A.InterpolationMethod = 'Flat';
AFlaIM = evalt(A,TGrid);

% Restore Warning Status
warning('on','ltvtools:evalt:extrapolate');

%% Plot
figure(1);clf;grid on;hold on;box on;
plot(Tm.Time,Tm,'-k','LineWidth',1);
plot(AFlaIM.Time,AFlaIM,'.-c','LineWidth',2.5);
plot(ANeaIM.Time,ANeaIM,'.-g','LineWidth',2.5);
plot(ALinIM.Time,ALinIM,':.r','LineWidth',2.5);
plot(AsplIM.Time,AsplIM,'--b','LineWidth',2.5);
plot(A.Time,A,'o','Color','m','MarkerSize',8,'MarkerFaceColor','m');
legend('True','Flat','Nearest','Linear','Spline','Samples','FontSize',14);
ylabel('Data','FontSize',16);xlabel('Time (s)','FontSize',16);
title('Interpolation Methods','FontSize',16);xlim([0,2*pi]);