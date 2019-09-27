function Beta = getTVUncLevel(T0,Ts,Tf,UL1,UL2)
Beta{1} = LOCALTVBound(1,T0,Ts,Tf,UL1,UL2);
Beta{2} = LOCALTVBound(1,T0,Ts,Tf,UL2,UL1);
Beta{3} = LOCALTVBound(2,T0,Ts,Tf,UL1,UL2);
Beta{4} = LOCALTVBound(2,T0,Ts,Tf,UL2,UL1);
Beta{5} = LOCALTVBound(3,T0,Ts,Tf,UL1,UL2);
Beta{6} = LOCALTVBound(3,T0,Ts,Tf,UL2,UL1);
Beta{7} = LOCALTVBound(4,T0,Ts,Tf,UL1,UL2);
Beta{8} = LOCALTVBound(4,T0,Ts,Tf,UL2,UL1);
end

function Beta = LOCALTVBound(type,T0,Ts,Tf,ULevel1,ULevel2)
%% Input Processing
narginchk(4,6);
nin = nargin;
if nin == 0
    type = 1;
end
if nin == 4
    % Uncertainty Levels
    ULevel1 = 0.1;
    ULevel2 = 0.9;
end

%% Load Example Data
% load('FHSynEx.mat','T0','Tf','Ts');

%% Time-Varying SISO Uncertainty Construction
switch type
    case 1
        %% Square Pulse (Rapid Transition)
        T1 = round(T0+(Tf-T0)*0.1);
        Tgrid1 = T0:Ts:T1;
        N1 = length(Tgrid1);
        Tgrid2 = T1+Ts:Ts:Tf;
        N2 = length(Tgrid2);
        Beta = tvmat([ULevel1*ones(1,N1) ULevel2*ones(1,N2)],[Tgrid1 Tgrid2]);
    case 2
        %% Square Pulse (Rapid Transition)
        T1 = round(T0+(Tf-T0)*0.9);
        Tgrid1 = T0:Ts:T1;
        N1 = length(Tgrid1);
        Tgrid2 = T1+Ts:Ts:Tf;
        N2 = length(Tgrid2);
        Beta = tvmat([ULevel1*ones(1,N1) ULevel2*ones(1,N2)],[Tgrid1 Tgrid2]);
    case 3
        %% Smooth Transition
        T1 = round(T0+(Tf-T0)*0.45);
        T2 = round(Tf-(Tf-T0)*0.45);
        ULData = [ULevel1 ULevel1 ULevel1 ULevel2 ULevel2 ULevel2];
        Time = [T1 T1+0.01 T1+0.02 T2-0.02 T2-0.01 T2];
        UL0 = evalt(tvmat(ULevel1),T0:Ts:T1);
        UL1 = evalt(tvmat(ULData,Time,'Spline'),T1:0.01:T2);
        UL1.InterpolationMethod = 'Linear';
        UL2 = evalt(tvmat(ULevel2),T2:Ts:Tf);
        Beta = tvmerge(UL0,UL1,UL2);
    case 4
        %% Gaussian
        T1 = round(T0+(Tf-T0)*0.36);
        T2 = round(Tf-(Tf-T0)*0.36);
        ULData = [ULevel1*ones(1,length(T0:Ts:T1))...
            ULevel1 ULevel1...
            ULevel2...
            ULevel1 ULevel1...
            ULevel1*ones(1,length(T2:Ts:Tf))];
        Time = [T0:Ts:T1 T1+0.01 T1+0.02 (T1+T2)/2 T2-0.02 T2-0.01 T2:Ts:Tf];
        Beta = evalt(tvmat(ULData,Time,'Spline'),T0:Ts:Tf);
        Beta.InterpolationMethod = 'Linear';
end
end