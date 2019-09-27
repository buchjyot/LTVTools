function [Y,I] = min(X,varargin)

if nargin>=2 && isa(varargin{1},'tvmat')
    Y = varargin{1};
    [Y,I] = tvbinop(@min,X,Y,varargin{2:end});
else
    [Y,I] = tvunopfl(@min,X,varargin{:});
end