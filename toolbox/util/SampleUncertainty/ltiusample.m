function out = ltiusample(NormBnd,NSamples)
% Given the normbound, this function samples stable first order LTI
% uncertainties within the frequency range 0.1 to 100 rad/s

% Input Checking
narginchk(1,2);
nin = nargin;
switch nin
    case 1
        NSamples = 100;
end

% The output is NSamples x 1 cell array of dynamic systems.
out = cell(NSamples,1);

% Range of frequencies
w1 = 1;
w2 = 100;
w = w1 + (w2-w1).*rand(NSamples,1);
sgn = sign(-1+2*rand(NSamples,1));

% LTI Uncertainties
for i = 1:NSamples
    out{i} = sgn(i)*NormBnd*tf([1 -w(i)],[1 w(i)]);
end
end