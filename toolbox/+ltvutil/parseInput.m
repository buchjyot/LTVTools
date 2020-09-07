function varargout = parseInput(Vin,Vclass)
%% ParseInput
% This helper function sorts the input by class
%
% Example:
% [a,b,c] = ltvutil.parseInput(varargin,{'udyn','double','tvwcOptions'});

nC = length(Vclass);
varargout = cell(nC,1);
for i = 1:nC
    idx = cellfun(@(x) isa(x,Vclass{i}), Vin);
    if sum(idx) == 1
        varargout{i} = Vin{idx};
    end
end
end