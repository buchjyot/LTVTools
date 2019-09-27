function H = plot(varargin)

% XXX Does not handle constant TVMATs

% Replace each TVMAT input by its (reshaped) data
nin = nargin;
for i=1:nin
    vi = varargin{i};
    if isa(vi,'tvmat')
        % Shift TVMAT Data so that first dimension corresponds to Time
        viData = vi.Data;
        nd = ndims(viData);
        viData = shiftdim(viData,nd-1);
        varargin{i} = viData;
    end
end

if nargout>0
    H = plot(varargin{:});
else
    plot(varargin{:});
end

