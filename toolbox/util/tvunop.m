function varargout = tvunop(opfh,A,varargin)
%% TVUNOPFL Unary operations in a for-loop for a TVMAT
%
% TVUNOP is a utility function that handles unary operations (e.g. transpose,
% ctranspose, etc.) on TVMATs. The user should not call TVUNOP directly.
% Instead each unary operation has a front end function that the user
% calls (e.g. transpose). The front end function then call calls
% TVUNOP  to do the required calculations/operations.

%% Initialize Output
nout = nargout;
varargout = cell(1,nout);

%% Perform unary operation
AData = A.Data;
if A.isTimeInvariant
    [varargout{1:nout}] = opfh(AData,varargin{:});
    for k=1:nout
        tmp = varargout{k};
        if isa(tmp,'double') || isa(tmp,'logical')
            varargout{k} = tvmat(tmp);
        elseif isa(tmp,'ss') || isa(tmp,'tf') || isa(tmp,'zpk')
            varargout{k} = tvss(tmp);
        elseif isa(tmp,'umat')
            varargout{k} = tvumat(tmp);
        elseif isa(tmp,'uss')
            varargout{k} = tvuss(tmp);
        end
    end
else
    szA = size(AData);
    nad = numel(szA)-3;
    id = repmat({':'},1,nad);
    Nt = numel(A.Time);
    out = cell(Nt,nout);
    for i=1:Nt
        [out{i,:}] = opfh( AData(:,:,id{:},i), varargin{:} );
    end
    for k=1:nout
        tmp = cat( ndims(out{1,k})+1 ,out{:,k} );
        if isa(tmp,'double') || isa(tmp,'logical')
            if isequal(A.Ts,0)
                varargout{k} = tvmat(tmp,A.Time,A.InterpolationMethod);
            else
                varargout{k} = tvmat(tmp,A.Time,A.Ts);
            end
        elseif isa(tmp,'ss') || isa(tmp,'tf') || isa(tmp,'zpk')
            varargout{k} = tvss(tmp,A.Time,A.InterpolationMethod);
        elseif isa(tmp,'umat')
            varargout{k} = tvumat(tmp,A.Time,A.InterpolationMethod);
        elseif isa(tmp,'uss')
            varargout{k} = tvuss(tmp,A.Time,A.InterpolationMethod);
        end
    end
end
