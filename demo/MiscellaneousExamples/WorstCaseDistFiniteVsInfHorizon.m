%% Horizon
T0 = 0;
Tf = 500;

%% System
Gtf = tf(10,[1 0.1 100]);
G = tvss(Gtf);

%% Compute hinfnorm and frequency
[g1,w] = hinfnorm(Gtf);
figure; bodemag(Gtf);

%% Compute finite horizon induced L2 gain
[g2,d,info] = tvnorm(evalt(G,[T0,Tf]));

%% Plotting
% Note the following:
% 1) The finite-horizon disturbance is sinusoidal
% 2) The disturbance frequency is approximately same as w
figure;
tvplot(d);