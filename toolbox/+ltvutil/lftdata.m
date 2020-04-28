function [M,Del,Blkstruct,NormUnc] = lftdata(obj)
%% LFTDATA Engine

if isempty(obj)
    [M,Del,Blkstruct,NormUnc] = lftdata(obj.Data);
    M = LOCALConvert2tv(M);Del = LOCALConvert2tv(Del);
    return;
end

if obj.isTimeInvariant
    [M,Del,Blkstruct,NormUnc] = lftdata(obj.Data);
    M = LOCALConvert2tv(M); Del = LOCALConvert2tv(Del);
else
    Nt = length(obj.Time);
    % XXX Preallocate Memory?
    for i = 1:Nt
        [MData(:,:,i),DelData(:,:,i),Blkstruct,NormUnc] = lftdata(obj.Data(:,:,:,i));
    end
    M = LOCALConvert2tv(MData,obj.Time,obj.InterpolationMethod);
    Del = LOCALConvert2tv(DelData,obj.Time,obj.InterpolationMethod);
end
end

function [M] = LOCALConvert2tv(Data,varargin)
switch class(Data)
    case 'ss'
        M = tvss(Data,varargin{:});
    case 'double'
        M = tvmat(Data,varargin{:});
    case 'umat'
        M = tvumat(Data,varargin{:});
    case 'uss'
        M = tvuss(Data,varargin{:});
end
end