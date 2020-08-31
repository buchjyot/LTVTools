%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2019                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implement a finite time horizon LQR controller for the Sit-to-Stand 
% movement of a minimally actuated Powered Lower Limb Orthosis at the hips. 
% The complete algorithm is published in:
% [1] O. Narvaez-Aroche, A. Packard, and M. Arcak, “Finite time robust 
% control of the Sit-to-Stand movement for powered lower limb orthoses,” 
% American Control Conference, June 2018, Milwaukee, WI, USA.
% http://dx.doi.org/10.23919/ACC.2018.8431465
%
% The reachability analysis under parameter uncertainty of the system is
% performed in: 
% [2] O.Narvaez-Aroche, P.-J.Meyer, M. Arcak, and A. Packard, “Reachability 
% analysis for robustness evaluation of the Sit-to-Stand movement for 
% powered lower limb orthoses,” ASME 2018 Dynamic Systems & Control 
% Conference, October 2018, Atlanta, GA, USA.
% http://dx.doi.org/10.23919/10.1115/DSCC2018-9066
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set a clean workspace
clear
clc
close all

%% Simulation parameters. 

% Total duration of ascension phase of the STS movement. 
tf = 3.5; % [s].

% Sitting position on the z-space for STS 1. 
zi = [-90*pi/180; 0.3099; 0.6678];

% Standing position on the z-space. 
zf = [-5*pi/180; 0; 0.9735];

% Nominal values for parameters of the system.
pnom = [9.68; 12.59; 44.57; 1.165; 0.519; 2.557; 0.533; 0.406; 0.52; 0.533/2; 0.406/2; 0.52/2];

% A different value for the parameters of the system.
punc = pnom + 0.01;

% Fixed initial state.
x0 = [90; -90; 90; 0; 0; 0]*pi/180;

%% Bounds for parameters uncertainty in reference [2]

% Number of parameters.
np = numel(pnom);

% Total nominal mass of links.
mT = sum(pnom(1:3));

% Total mass of links variation factor.
mvf = 0.05;

% Arrays for parameter bounds.
pmax = zeros(1,np);
pmin = zeros(1,np);

% Apply uncertainties on mass of links [kg].
pmax(1:3) = [(1+mvf)*0.1448*mT, (1+mvf)*0.1884*mT, (1+mvf)*0.6668*mT];
pmin(1:3) = [(1-mvf)*0.1448*mT, (1-mvf)*0.1884*mT, (1-mvf)*0.6668*mT];

% Apply uncertainties on length of links [m].
pmax(7:9) = pnom(7:9) + 0.01;
pmin(7:9) = pnom(7:9) - 0.01;

% Apply uncertainties on the position of the CoM of each link [m].
pmax(10:12) = [0.567*pnom(7), 0.567*pnom(8), (pnom(9)/2)+0.0154];
pmin(10:12) = [0.433*pnom(7), 0.433*pnom(8), (pnom(9)/2)-0.0154];

% Apply uncertainties on moments of inertia of links [kg.m^2].
pmax(6) = (1/12)*pmax(3)*((0.81)^2+(0.15)^2);
pmin(6) = (1/12)*pmin(3)*((0.81)^2+(0.15)^2);
pmax(4:5) = [(1+(pmax(6)-pnom(6))/pnom(6))*pnom(4), (1+(pmax(6)-pnom(6))/pnom(6))*pnom(5)];
pmin(4:5) = [(1-(pnom(6)-pmin(6))/pnom(6))*pnom(4), (1-(pnom(6)-pmin(6))/pnom(6))*pnom(5)];

%% Integrate systems dynamics under LTV LQR Control over one time step.

% Load array with LQR gains and corresponding time scale. 
load('LQRGainSTS1.mat')

% Bounds for states during ode45 simulation of system under TV LQR control.
xmin = [70;-130;-10;-25;-20;-80]*pi/180;
xmax = [130;10;150;15;70;30]*pi/180;

% Set options for ode45 events function.
options = odeset('Events',@(t,x)TVLQRSTS3LinkEvents(t,x,xmin,xmax));

% Solve ODE.
fprintf('\nPerforming Time-Varying LQR control for STS 1 with nominal parameter...\n')
tic
[Tnom,Xnom] = ode45(@(t,x) LQRThreeLink(t,tf,x,KLQRSTS1,tgrid,zi,zf,pnom,pnom),[0 tf],x0,options);
toc

% Solve ODE for the lower bound of the parameter uncertainty. 
fprintf('\nPerforming Time-Varying LQR control for STS 1 with lower parameter bound...\n')
tic
[Tmin,Xmin] = ode45(@(t,x) LQRThreeLink(t,tf,x,KLQRSTS1,tgrid,zi,zf,pnom,pmin),[0 tf],x0,options);
toc

% Solve ODE for the upper bound of the parameter uncertainty. 
fprintf('\nPerforming Time-Varying LQR control for STS 1 with upper parameter bound...\n')
tic
[Tmax,Xmax] = ode45(@(t,x) LQRThreeLink(t,tf,x,KLQRSTS1,tgrid,zi,zf,pnom,pmax),[0 tf],x0,options);
toc

%% Obtain output trajectories.

% Reference.
Ynom = x2yThreeLink(Xnom',pnom)';

% For lower bound of the parameter.
Ymin = x2yThreeLink(Xmin',pmin)';

% For upper bound of the parameter.
Ymax = x2yThreeLink(Xmax',pmax)';

%% Obtain input trajectories.

% Reference.
Unom = LQRInputThreeLink(Tnom,Xnom,KLQRSTS1,tgrid,zi,zf,pnom);

% For lower bound of the parameter.
Umin = LQRInputThreeLink(Tmin,Xmin,KLQRSTS1,tgrid,zi,zf,pnom);

% For upper bound of the parameter.
Umax = LQRInputThreeLink(Tmax,Xmax,KLQRSTS1,tgrid,zi,zf,pnom);

%% Plot simulations.

% State.
figure()
for i=1:numel(x0)
   subplot(2,3,i)
   plot(Tnom,Xnom(:,i)*180/pi,'r-',Tmin,Xmin(:,i)*180/pi,'b-',Tmax,Xmax(:,i)*180/pi,'g-','LineWidth',2)
   grid()
   xlabel('$t \, [s]$','Interpreter','Latex')
   if i<=3
       ylabel(['$\theta_{',num2str(i),'}(t) \, [^\circ]$'],'Interpreter','Latex')
   else
       ylabel(['$\dot{\theta}_{',num2str(i-3),'}(t) \, [^\circ/s]$'],'Interpreter','Latex')
   end
   legend({'$p = \hat{p}$','$ p = \underline{p}$','$ p = \overline{p}$'},'Interpreter','Latex','Location','Best')
   set(gca,'FontSize',14)
   xlim([Tnom(1),Tnom(end)])
end

% Output.
figure()
for i = 1:6
    subplot(2,3,i)
    if i<5
        plot(Tnom,Ynom(:,i),'r-',Tmin,Ymin(:,i),'b-',Tmax,Ymax(:,i),'g-','LineWidth',2)
        xlim([tgrid(1), tgrid(end)]);
        xlabel('$t \, [s]$','Interpreter','Latex');
    elseif i==5
        plot(Ynom(:,1),Ynom(:,2),'r-',Ymin(:,1),Ymin(:,2),'b-',Ymax(:,1),Ymax(:,2),'g-','LineWidth',2)
    else
        plot(Ynom(:,3),Ynom(:,4),'r-',Ymin(:,3),Ymin(:,4),'b-',Ymax(:,3),Ymax(:,4),'g-','LineWidth',2)
    end
    grid()
    if i==1
        ylabel('$ x_{CoM}(t) \, [m] $','Interpreter','Latex');
    elseif i==2
        ylabel('$ y_{CoM}(t) \, [m] $','Interpreter','Latex');
    elseif i==3
        ylabel('$ \dot{x}_{CoM}(t) \, [m/s] $','Interpreter','Latex');
    elseif i==4
        ylabel('$ \dot{y}_{CoM}(t) \, [m/s] $','Interpreter','Latex');
    elseif i==5
        xlabel('$ x_{CoM}(t) \, [m] $','Interpreter','Latex');
        ylabel('$ y_{CoM}(t) \, [m] $','Interpreter','Latex');
    else
        xlabel('$\dot{x}_{CoM}(t) \, [m/s]$','Interpreter','Latex');
        ylabel('$\dot{y}_{CoM}(t) \, [m/s]$','Interpreter','Latex');
    end
    legend({'$p = \hat{p}$','$ p = \underline{p}$','$ p = \overline{p}$'},'Interpreter','Latex','Location','Best')
    set(gca,'Fontsize',14)
end

% Input.
figure()
for i = 1:4
    subplot(2,2,i);
    plot(Tnom,Unom(:,i),'r-',Tmin,Umin(:,i),'b-',Tmax,Umax(:,i),'g-','LineWidth',2)
    xlim([tgrid(1), tgrid(end)]);
    grid()
    % Add labels.
    xlabel('$t \, [s]$','Interpreter','Latex');
    if i==1
        ylabel('$\tau_{h}(t) \, [N.m]$','Interpreter','Latex')
    elseif i==2
        ylabel('$\tau_{s}(t)\,[N.m]$','Interpreter','Latex')
    elseif i==3
        ylabel('$F_{x}(t)\,[N]$','Interpreter','Latex')
    else
        ylabel('$F_{y}(t)\,[N]$','Interpreter','Latex')
    end
    legend({'$p = \hat{p}$','$ p = \underline{p}$','$ p = \overline{p}$'},'Interpreter','Latex','Location','Best')
    set(gca,'Fontsize',14)
end