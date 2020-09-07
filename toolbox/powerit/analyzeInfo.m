function fh = analyzeInfo(info)
%% ANALYZE RESULTS
% This function analyzes the results based on the input info structure

% Read Info Structure
Niter = info.TotalIter;
if isfield(info,'AllIter')
    AllIter = info.AllIter;
else
    fh = [];
    return
end
iterVec = 1:Niter;
FwdPerf = [AllIter.FwdPerf];
AdjPerf = [AllIter.AdjPerf];
Nwc = info.WCIter;

% Plot forward performance vs iteration counts
fh(1) = figure;
plot(iterVec(:),FwdPerf(:),'s--','LineWidth',1);
xlabel('Iteration Count','FontSize',14);
ylabel('Objective','FontSize',14);
grid on;box on;
set(gca,'GridLineStyle','--');
xl = xlim;xlim([1,xl(2)]);

% Plot performance vs iteration counts
fh(2) = figure;
plot(iterVec(:),FwdPerf(:),'s--b',iterVec(:),AdjPerf(:),'o--r','LineWidth',1);
xlabel('Iteration Count','FontSize',14);
ylabel('Objective','FontSize',14);
grid on;box on;
set(gca,'GridLineStyle','--');
legend('Forward Performance','Adjoint Performance','Location','southeast','FontSize',14);
xl = xlim;xlim([1,xl(2)]);

% Check if iterations were monotincally improving
CheckMonotone = false;
if CheckMonotone 
    if all(diff(FwdPerf(:)) >= -1e-2)
        monotonicString = 'Monotonic';
    else
        monotonicString = 'Non-monotonic';
    end
    title(sprintf('Iteration Progress: %s',monotonicString),'FontSize',14);
end

% Plot sequence of input signals
fh(3) = figure;
for k = 1:Niter
    if k~=Nwc
        thisU = AllIter(k).U;
        tvplot(thisU,'Color','g','LineWidth',1);
        hold on;
    end
end
tvplot(AllIter(Nwc).U,'Color','b','LineWidth',1.5);
xlabel('Time (s)','FontSize',14);
ylabel('Input','FontSize',14);
title('Sequence of Inputs','FontSize',14);
grid on;box on;
set(gca,'GridLineStyle','--')

% Plot worst-case input signals only
fh(4) = figure;
tvplot(AllIter(Nwc).U,'Color','b','LineWidth',1.5);
xlabel('Time (s)','FontSize',14);
ylabel('Input','FontSize',14);
title('Worst-Case Input','FontSize',14);
grid on;box on;
set(gca,'GridLineStyle','--')
end