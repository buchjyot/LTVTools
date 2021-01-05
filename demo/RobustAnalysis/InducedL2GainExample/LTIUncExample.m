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

% Form uncertain system
Delta = ultidyn('Delta',[1 1]);
Gunc = lft(Delta,G);

%% MU/SSV: Infinite-Horizon Robust Stability and Worst-Case Gain
[StabMarg,WCU1] = robstab(Gunc);
[WCGain,WCU2] = wcgain(Gunc);
gWCG = WCGain.UpperBound;

fprintf('\n Stab Margin (LB) = %4.3f at w = %4.3f',...
    StabMarg.LowerBound,StabMarg.CriticalFrequency);
fprintf('\n WC Gain (UB) = %4.3f at w = %4.3f',...
    gWCG,WCGain.CriticalFrequency);

%% IQC: Infinite-Horizon Worst-Case Gain
v = 1;
p = -10;  % -1, -10
[gih,X11ih,Pih] = IHL2toL2lmi(G,v,p);
fprintf('\n IH IQC Gain Bound = %4.3f\n',gih);

%% Robust L2 to L2 Gain on Finite Horizon
% Iterate between finite horizon LMI and RDE Solutions

% Perfom Analysis over a range of horizons
Tall = [1 2 5 10:10:50 100]; %200];
NT = numel(Tall);
AllResults = struct('T',{},'wcinfo',{},'gfinal',{},'Niter',{});
n = size(Pih,1);

% Analysis
Nall = zeros(NT,1);
RobTime = zeros(NT,1);
DeltaIQC = udyn('Delta',[1 1],'UserData',[0,p,v]);
tvwcopt = tvwcOptions('Nlmi',20,'Nsp',10,'Display','on');
T0 = 0;
tstart = tic;
for k=1:NT
    % Finite Horizon
    T = Tall(k);
    fprintf('\n Analysis for T=%4.3f \t (%d of %d)',T,k,NT);
    
    % Robust L2 to L2 gain
    % Deprecated Syntax: [gfinal,wcinfo] = tvrobL2toL2(tvss(G),v,p,T,tlmi,tSp);
    Gtv = tvss(G,[T0 T]);
    [gfinal,wcinfo] = tvwcgain(Gtv,DeltaIQC,0,tvwcopt);
    Niter = numel(wcinfo.AllIter);
    Nall(k) = Niter;
    RobTime(k) = wcinfo.TotalTime;
    
    % Store final Results
    AllResults(k) = struct('T',T,'wcinfo',wcinfo,'gfinal',gfinal,'Niter',Niter);
end
tstop = toc(tstart);
fprintf(newline);

%% Nominal Analysis
Gnom = lft(0,G);
[Ginf,winf] = norm(Gnom,inf);
fprintf('\n (IH) Nominal Gain = %4.3f at w = %4.3f',Ginf,winf);
gT = zeros(NT,1);
NomTime = zeros(NT,1);
for k=1:NT
    % Finite Horizon
    T = Tall(k);
    fprintf('\n Analysis for T=%4.3f \t (%d of %d)',T,k,NT);
    
    tic;
    gg = tvnorm(evalt(tvss(Gnom),[T0 T]));
    gT(k) = gg(2);
    NomTime(k) = toc;
    fprintf('\n FH Nom Gain = %4.3f',gT(k));
end
fprintf('\n');

%% Plot Results
gAll = [AllResults.gfinal];

if numel(Tall)>1
    figure;
    ph1a = plot(Tall,gAll,'b',Tall,gAll,'bx',Tall([1 end]),gih*[1 1],'r--');
    %hold on;
    %ph1b = plot(Tall,gT,'r',Tall,gT,'rx',Tall(end),Ginf,'ro');
    %hold off
    xlim([0 100])
    xlabel('Horizon Time T, sec')
    ylabel('L_2 Gain');
    grid on;
    
    figure
    plot(Tall,RobTime./Nall*2,'b',Tall,RobTime./Nall*2,'bx');
    hold on;
    plot(Tall,NomTime,'r',Tall,NomTime,'rx');
    hold off
    xlabel('Horizon Time T, sec')
    ylabel('Comp Time (sec)');
    grid on;
end

return
%% Simple computational study
% This code performs a simple comparison of the computation required
% for different approaches to assess the finite-horizon robustness.
% Results are for the horizon T=20 using the data above.
%
% A) Proposed algorithm - The results are printed below using the settings
%    Nsp = 10 and Nlmi= 20.
%      Iteration # = 1
%      LMI Gain Bound = 1.3062,	 RDE Gain Bound = 1.2197,
%      Iteration # = 2
%      LMI Gain Bound = 1.2190,	 RDE Gain Bound = 1.2168,
%      Final Results: Robust gain = 1.2168  Total Comp Time = 39.1213
%
% B) "Approximate" LMI - Solve the problem using spline bases with
%    Nsp time points and an LMI time grid with 2*Nsp time points.
%    Results are given below for several values of Nsp
%       Nsp = 10  gFinal = 1.3062   tStop = 6.0438
%       Nsp = 20  gFinal = 1.2542   tStop = 68.6675
%       Nsp = 30  gFinal = 1.2399   tStop = 401.0048
%    Note that the results for Nsp=10 correspond to the first LMI
%    calculation of the proposed algorithm (A).
%
% C) RDE - We implemented a simple cutting plane algorithm for a single
%    IQC variable. We did not implement a more general RDE-based cutting
%    plane algorithm capable of handling many IQC variables. However, we
%    evaluated the computational cost of using the RDE.
%    Specifically, the computation to evaluate the cost using RDE +
%    bisection (with AbsTol = 1e-4 and RelTol = 1e-3) was 11.1578 secs.
%    This cost evaluation required 12 calls to CDRE.  Thus each RDE
%    took ~0.93sec on this horizon.  The total comp time for the proposed
%    algorithm A corresponds to roughly 3.5 function evals (=39.12/11.16)
%    using RDE + bisection.  It also corresponds to roughly 42 RDE
%    evals (=39.12/0.93).  It would be difficult for an RDE-based approach
%    to be computationally competitive if there is more than one IQC var.


% Choose time horizon
T = 20;

% Use Nsp points for spline and 2*Nsp for lmi
Nsp = 20;
tlmi = linspace(0,T,2*Nsp);
tSp = linspace(0,T,Nsp);
Ps = [];
for i=1:Nsp
    ei = zeros(Nsp,1);
    ei(i) = 1;
    Ps = [Ps; tvmat(ei,tSp,'Spline')];
end

% Finite Horizon: LMI Condition
G1 = evalt(tvss(G),tlmi);
Ps1 = evalt(Ps,tlmi);
Psdot1 = tvdiff(Ps,tlmi);
Pm = tvmat; Pmdot = Pm;

tStart = tic;
[gFinal,X11]= FHL2toL2lmi(G1,v,p,Ps1,Psdot1,Pm,Pmdot);
gFinal
tStop = toc(tStart)

% Finite Horizon: RDE
tStart2 = tic;
Nrde = 5;
for i=1:Nrde
    [gbnds,RDEinfo] = FHL2toL2rde(tvss(G),v,p,T,X11);
end
tAvg = toc(tStart2)/Nrde
RDEinfo.RDEcnt