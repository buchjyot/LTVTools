function samples = uniformsampling(Nom,Del,N)
% Inputs:
%
% Nom := NominalValue
% Del := Range of PercentageUncertainty
% N   := Number of Samples

% e.g. [-10, 10] Means -10% to +10%
Del1 = Del(1);
Del2 = Del(2);

LB = Nom + (Nom*Del1/100);
UB = Nom + (Nom*Del2/100);

samples = LOCALIntervalSampling(LB,UB,N);
end

function out = LOCALIntervalSampling(a,b,N)
out = a + (b-a).*rand(N,1);
end