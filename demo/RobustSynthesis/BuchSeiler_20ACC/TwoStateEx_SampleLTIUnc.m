%% Load Data
load('TwoStateEx.mat');

%% LTI Uncertainty
gain = [linspace(-beta,-beta+0.1,50) linspace(beta-0.1,beta,50)];
N = length(gain);
UncSample = cell(N+2,1);
for i = 1:N
    b = abs(randn)*10;
    UncSample{i} = gain(i)*tf([1 -b],[1 b]);
end

%% Constant Uncertainty
UncSample{N+1} = beta;
UncSample{N+2} = -beta;

%% Save Data
save(mfilename,'UncSample');