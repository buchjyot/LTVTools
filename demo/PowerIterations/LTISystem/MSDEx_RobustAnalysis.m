%% Mass-Spring-Damper System
% m := mass
% b := damping coefficient
% k := spring constant
msdEx1 = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
P  = msdEx1(1,2*0.7,1);

% Transform plant to have all the states (Nx) as outputs
Nx = order(P);
P = ss(P.A,P.B,eye(Nx),0);
[Ne,Nd] = size(P);

%% Horizon
T0   = 0;
Tall = 3;%[0.1,1,2,3,4,5,7,10,20,50,70];
NT   = length(Tall);

%% Analysis Interconnection
% Input Scalling
WdScl = 0.1;

% Output Scalling
WexScl = 1;

% Uncertainty Scalling
Nv = 1;
Nw = 1;
DelNorm = 0.6;
DelInScl = sqrt(DelNorm);
DelOutScl = sqrt(DelNorm);

% Promote to state-space objects
Wd = ss(WdScl*eye(Nd));
We = ss(WexScl*eye(Ne));
WDelIn = ss(DelInScl*eye(Nv));
WDelOut = ss(DelOutScl*eye(Nw));

% System Interconnection
systemnames = 'P Wd We WDelIn WDelOut'; %#ok<*NASGU>
inputvar = '[w; d]';
outputvar = '[WDelIn; We]';
input_to_P = '[WDelOut+Wd]';
input_to_Wd = '[d]';
input_to_We = '[P]';
input_to_WDelIn = '[Wd]';
input_to_WDelOut = '[w]';
cleanupsysic = 'yes';
Gunc = sysic;

% Euclidean Penalties
NE = 0;

%% Uncertainty
DeltaIQC = udyn('Del',[1 1],'UserData',[0,-10,1]);
Delta = ultidyn('Del',[1 1],'Bound',1);

%% Options
OdeSolver = 'ode23s';

pOpt  = poweritOptions('Display','on','OdeSolver',OdeSolver,'StoreAllIter',true);
pSpec = poweritSignalSpec('Nv',1,'Nw',1,'NE',NE,'InitialInput','randn');
wcOpt = tvwcOptions('Display','on');
wcOpt.RDEOptions.OdeSolver = OdeSolver;

% seed
rng(0);

%% Infinite Horizon Gains
if isequal(NE,0)
    [wcgIH,wcuIH,wcInfoIH] = wcgain(lft(Delta,Gunc));
    fprintf(' Infinite Horizon wcgLB: %.4f, wcgUB: %.4f\n\n',wcgIH.LowerBound,wcgIH.UpperBound);
end

%% Worst-Case Gain Upper Bounds
ComputeWcGainUpperBounds = false;
if ComputeWcGainUpperBounds
    wcgUB = zeros(NT,1);
    wcInfo = cell(NT,1);
    for i = 1:NT
        [wcgUB(i),wcInfo{i}] = tvwcgain(tvss(Gunc,[T0,Tall(i)]),DeltaIQC,NE,wcOpt);
        fprintf(newline);
    end
end

%% Power Iterations Lower Bound
ComputeWcGainLowerBounds = true;
if ComputeWcGainLowerBounds
    wcgLB = zeros(NT,1);
    pInfo = cell(NT,1);
    wcSig = cell(NT,1);
    for i = 1:NT
        [wcgLB(i),wcSig{i},pInfo{i}] = powerit(Gunc,[T0,Tall(i)],pSpec,pOpt);
        fprintf(newline);
    end
end

%% Save Data
fprintf(' Bounds:\n');
disp([wcgLB wcgUB]);
save(mfilename,'Tall','wcgIH','wcuIH','wcInfoIH','wcgUB','wcInfo','wcgLB','wcSig','pInfo','Gunc','NE');

%% Plot Upper and Lower Bounds
if length(Tall) > 1
    figure(1);clf;
    plot(Tall,wcgLB,'-^b',Tall,wcgUB,'-.vr','LineWidth',2);
    grid on;box on;
    xlabel('Horizon (T) (sec.)');
    ylabel('Worst-case Gain');
    legend('Power Iteration Lower Bound','IQC Upper Bound','Location','southeast');
    xlim([0 25]);
end