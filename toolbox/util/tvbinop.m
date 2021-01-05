function C = tvbinop(opfh,A,B,varargin)
% TVBINOP Binary operations for a TVMAT
%
% TVBINOP is a utility function that handles binary operations (e.g. times,
% ldivide, lft, etc.) on TVMATs. The user should not call TVBINOP directly.
% Instead each binary operation has a front end function that the user
% calls (e.g. *,+,lft,horzcat). The front end function then call calls
% TVBINOP  to do the required calculations/operations.

%% Use switchyard to handle object interactions
[isTimeInvariant,A,B]=tvswitchyard(A,B);

%% Handle different TVMAT data types
if isTimeInvariant
    % Both A and B are constant
    CData = opfh( A.Data, B.Data, varargin{:} );
    if isa(CData,'double') || isa(CData,'logical')
        C = tvmat(CData);
    elseif isa(CData,'ss') || isa(CData,'tf') || isa(CData,'zpk')
        C = tvss(CData);
    elseif isa(CData,'umat')
        C = tvumat(CData);
    elseif isa(CData,'uss')
        C = tvuss(CData);
    end
else
    % Both A and B are Grids
    ATime = A.Time;
    Nt = numel(ATime);
    AIM = A.InterpolationMethod;
    ATs = A.Ts;
    
    % Grab TVMAT Data
    AData = A.Data;
    nad = ndims(AData)-3;
    aid = repmat({':'},1,nad);
    
    BData = B.Data;
    nbd = ndims(BData)-3;
    bid = repmat({':'},1,nbd);
    CData = cell(Nt,1);
    for i=1:Nt
        CData{i} = opfh( AData(:,:,aid{:},i), BData(:,:,bid{:},i), varargin{:} );
    end
    CData = cat( ndims(CData{1})+1, CData{:});
    if isa(CData,'double') || isa(CData,'logical')
        if isequal(ATs,0)
            C = tvmat(CData,ATime,AIM);
        else
            C = tvmat(CData,ATime,ATs);
        end
    elseif isa(CData,'ss') || isa(CData,'tf') || isa(CData,'zpk')
        C = tvss(CData,ATime,AIM);
    elseif isa(CData,'umat')
        C = tvumat(CData,ATime,AIM);
    elseif isa(CData,'uss')
        C = tvuss(CData,ATime,AIM);
    end
end
end