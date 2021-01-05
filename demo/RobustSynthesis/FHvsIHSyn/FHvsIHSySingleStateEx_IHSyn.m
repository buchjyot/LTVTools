%% Load Data
load('FHvsIHSySingleStateEx.mat');

%% Perform Infinite Horizon Robust Synthesis
Gu = lft(Del,Gunc);
Gscl = Gu;

gUpp = 2;
gLow = 0.01;
gTry = (gLow+gUpp)/2;
Tol = 1e-4;

while true
    % Scale performance channel until Robust Performance from DKsyn is 1
    Gscl = Gu*blkdiag(eye(NU-Nw-Ny)/gTry,eye(Ny));
    
    % Run DK Iterations
    [K,CL,GAM,DKINFO] = dksyn(Gscl,Ny,Nu,opt);
    gstr = dksynperf(CL);
    gRP  = gstr.UpperBound;
    
    % Bisect until RobustPerformance ~= 1
    fprintf('gUpp = %.4f, gLow = %.4f, gTry = %.4f, gRP = %.4f\n',gUpp,gLow,gTry,gRP);
    if abs(gRP-1) <= Tol
        break;
    end
    if gRP > 1
        gLow = gTry;
    else
        gUpp = gTry;
    end
    gTry = (gLow+gUpp)/2;
end
gWCGain = gTry;
fprintf('LTI Synthesis worst-case gain : %.4f\n',gWCGain);

%% Save
save(mfilename);