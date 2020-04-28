%% Problem Data
% LTI State-Space
G = ss([-1 0; 0 -2],[1;-1],[1,1],0);

% Transfer function
Nx = order(G);
Gtf = tf(G);

% Finite Horizon System
T0 = 0;
Ts = 0.01;
Tf = 5;
Opt = tvodeOptions('OdeSolver','ode23tb','OdeOptions',odeset('RelTol',1e-4,'AbsTol',1e-4));
Gfh = evalt(tvss(G),T0:Ts:Tf);

% Balanced Realization
[Gb,hSV,Ts,Tsi,info] = tvbalreal(Gfh,Opt);