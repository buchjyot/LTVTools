function verifyFH(varargin)
%% verifyFH - verifies that tvobject is finite horizon
nin = nargin;
for i = 1:nin
    if ~isfh(varargin{i})
        error('Horizon must be finite.');
    end
end
end