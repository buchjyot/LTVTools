%% TVKALMAN 2by2 Example
% Use tvkalman to design an LTV Kalman Filter for SISO 2by2 plant
%
% This example is taken from: 
% 
% Lewis, F.L. 1986. Optimal estimation: with an introduction to stochastic
% control theory (pp. 88-94). New York et al.: Wiley. Page 164, Example 3.4
% Damped Harmonic Oscillator

%% System Specification

% Specify Time Grid
T0 = 0;
Tf = 10;
TGrid = linspace(T0,Tf,1000);

omegan = 0.64;
alpha = 0.16;
b = 0;
A = [0 1;-omegan^2 -2*alpha];
B = [0; b];
G = [0; 1];
C = [0 1];
D = 0;
H = 0;

% Process Noise Covariance
Qn = 1;

% Measurement Noise Covariance
Rn = 1;

%% Finite Horizon Kalman Filter

% Time-Varying System for Kalman Filter Design
%      .
%      x = A(t)x(t) + B(t)u(t) + G(t)*w        {State equation}
%      y = C(t)x(t) + D(t)u(t) + H(t)*w + v    {Measurements}
A = evalt(tvmat(A),TGrid);
sys = tvss(A,[B G],C,[D H]);

P0 = zeros(2);
[Kfh,Lfh,Pfh] = tvkalman(sys,Qn,Rn,P0,[T0 Tf]);
[Ae,Be,Ce,De] = ssdata(Kfh);

% Plot covariances 
% Figure 3.8 in the book
figure(1);clf;
tvplot(Pfh,'LineWidth',2);
title('Error covariance terms for a harmonic oscillator.')
legend('P1(t)','P12(t)=P21(t)','P22(t)');
grid on;

% Plot kalman gains
figure(2);clf;
plot(Lfh,'LineWidth',2);
title('Kalman Gains');
legend('L1','L21');
grid on;