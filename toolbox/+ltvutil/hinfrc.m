function [r1,r2] = hinfrc(D12,D21,Ny,Nu)
%% HINFRC
% Utility function that checks the rank conditions for the hinfinity and h2
% synthesis routines.
eflag = false;

% Condition 1
rank1 = rank(D12'*D12);
r1 = unique(rank1.Data);
if ~isequal(r1,Nu) || ~isequal(numel(r1),1)
    eflag = true;
end

% Condition 2
if ~isequal(Ny,0)
    rank2 = rank(D21*D21');
    r2 = unique(rank2.Data);
    if ~isequal(r2,Ny) || ~isequal(numel(r2),1)
        eflag = true;
    end
end

% Throw error if rankcond is not satisfied
if eflag
    error('Rank conditions are not satisfied.');
end
end