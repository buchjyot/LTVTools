function out = tvrecop(opfh,varargin)
% TVRECOP N-ary operations for a TVMAT implented as recursive binary ops

% Check # of input arguments
narginchk(2, inf)

v1 = varargin{1};
if isa(v1,'double') || isa(v1,'logical')
    v1 = tvmat(v1);
elseif isa(v1,'ss') || isa(v1,'tf') || isa(v1,'zpk')
    v1 = tvss(v1);
elseif isa(v1,'umat')
    v1 = tvumat(v1);
elseif isa(v1,'uss')
    v1 = tvumss(v1);
end

if nargin==2
    out = v1;
else
    out = tvbinop(opfh,v1,varargin{2});
    if nargin>2
        out = opfh(out,varargin{3:end});
    end    
end

