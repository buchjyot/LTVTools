rng(0);

% Matrix
M = randn(3);
M = M/norm(M);

% Uncertainty
Del = ureal('Del',0);

% LFT
P = lft(Del,M);

% Worst-Case Norm
[wcn,wcu] = wcnorm(P)