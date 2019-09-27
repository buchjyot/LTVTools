%% Load Data
load('TwoStateEx.mat');

%% LTI Uncertainty
gain = [linspace(-beta,-beta+0.1,20) linspace(beta-0.1,beta,20)];
N = length(gain);
Delta = cell(N+2,1);
for i = 1:N
    b = abs(randn)*10;
    Delta{i} = gain(i)*tf([1 -b],[1 b]);
end

%% Constant Uncertainty
Delta{N+1} = beta;
Delta{N+2} = -beta;

%% Save Data
save(mfilename,'Delta');