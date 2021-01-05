function out = ltiusample(NormBnd,NSamples,MaxOrder)
% Given the normbound, this function samples stable, SISO, LTI uncertainty
% using RSS command, the order of Delta is anything between 1 and MaxOrder

% Input Checking
narginchk(0,3);
nin = nargin;
switch nin
    case 0
        % Default is single state, SISO, norm bounded by 1, uncertainties
        NormBnd = 1;
        NSamples = 100;
        MaxOrder = 1;
    case 1
        NSamples = 100;
        MaxOrder = 1;
    case 2
        MaxOrder = 1;
end

% The output is NSamples x 1 cell array of dynamic systems.
out = cell(NSamples,1);

% LTI Uncertainties
switch MaxOrder
    case 1
        % Range of frequencies
        w1 = 1;
        w2 = 100;
        w = w1 + (w2-w1).*rand(NSamples,1);
        sgn = sign(-1+2*rand(NSamples,1));
        
        % LTI Uncertainties
        for i = 1:NSamples
            out{i} = sgn(i)*NormBnd*tf([1 -w(i)],[1 w(i)]);
        end
        
    otherwise
        for i = 1:NSamples
            SYS = rss(randi(MaxOrder),1,1);
            hn = hinfnorm(SYS);
            while ~isstable(SYS) || ~isfinite(hn)
                SYS = rss(randi(MaxOrder),1,1);
                hn = hinfnorm(SYS);
            end
            out{i} = sqrt(NormBnd/hn)*SYS*sqrt(NormBnd/hn);
        end
end