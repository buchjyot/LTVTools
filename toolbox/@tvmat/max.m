function [Y,I] = max(X,varargin)

if nargin>=2 && isa(varargin{1},'tvmat')
    Y = varargin{1};
    [Y,I] = tvbinop(@max,X,Y,varargin{2:end});
else
    [Y,I] = tvunop(@max,X,varargin{:});
end