function samples = gaussiansampling(Mean,SD,N)
%% Sample from Gaussian Distribution
% Inputs:
%
% Mean := Nominal Mean Value
% SD   := Standard Deviation
% N    := Number of Samples

samples = Mean + sqrt(SD)*randn(N,1);
end