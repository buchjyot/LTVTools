function [A1,A2,A3] = tvdiff(A,T)
% TVDIFF Differentiate a TVMAT with respect to time.
%
% [A1,A2,A3] = TVDIFF(A,T) differentiates a TVMAT A at the times specified
% in the vector T. A1, A2, A3 are the first, second, and third derivatives 
% of A with respect to time. They are returned as TVMATs.  
%
% TVDIFF(A) is equivalent to TVDIFF(A,A.Time);
%
% TVDIFF requires A to have InterpolationMethod='Spline'. The derivative
% is computed by analytically using the cubic spline functions.

%% Check inputs
% XXX Currently requires IM = 'Spline'. 
% Implement approximate derivatives for other cases?
narginchk(1,2);
if ~isa(A,'tvmat')
    error('A must be a TVMAT');
end

AIM = A.InterpolationMethod;
if ~isequal(AIM,'Spline')
    error('TVDIFF requires Spline interpolation method');
end

if isdt(A)
    error('Discrete-Time TVMAT can not be differentiated.')
end

%% Handle single input syntax with recursive call
if nargin==1
    [A1,A2,A3] = tvdiff(A,A.Time);
    return   
end

%% Handle Constant Case
AData = A.Data;
Nd = ndims(AData);
Nt = numel(T);
if A.isTimeInvariant
    AdData = repmat( zeros(size(AData)) ,[ones(1,Nd) Nt]);
    A1 = tvmat(AdData,T);
    A2 = A1;
    A3 = A1;    
    return
end

%% Handle Varying Case
ATime = A.Time;

% Pre-process TVMAT Data
nad = ndims(AData)-3;
id = repmat({':'},1,nad);

% Currently SplineData is empty and it computed in every call
% to EVALT. A syntax to precompute Spline Data is:
%    A = getSplineData(A);
% This allows for more efficient, repeated calls to EVALT.
ASplineData = A.SplineData;
if isempty(ASplineData)
    ASplineData = getSplineData(AData,ATime);
end
        
% XXX The code below allows time data outside Time range but
% clips to nearest time point. Error out instead?
for i=1:Nt
    if T(i)<=ATime(1)
        k=1;
        dt = 0;
    elseif T(i)>=ATime(end)
        k=numel(ATime)-1;
        dt = ATime(end)-ATime(end-1);
    else
        k = find( T(i)>=ATime(1:end-1) & T(i)<ATime(2:end) );
        k = k(1);
        dt = T(i)-ATime(k);
    end
    aidx = 4*(k-1)+(1:4);
    if dt==0
        A1i = ASplineData(aidx(2),:);
        A2i = 2*ASplineData(aidx(3),:);
        A3i = 6*ASplineData(aidx(4),:);        
    else
        A1i = [0 1 2*dt 3*dt^2]*ASplineData(aidx,:);
        A2i = [0 0 2 6*dt]*ASplineData(aidx,:);
        A3i = [0 0 0 6]*ASplineData(aidx,:);
    end
    
    if i==1
        A1Data = zeros( [numel(A1i) Nt] );
        A2Data = zeros( [numel(A1i) Nt] );
        A3Data = zeros( [numel(A1i) Nt] );
    end
    A1Data(:,i) = A1i(:);
    A2Data(:,i) = A2i(:);
    A3Data(:,i) = A3i(:);
end
A1Data = reshape(A1Data,[size(AData(:,:,id{:},1)) Nt]);
A2Data = reshape(A2Data,[size(AData(:,:,id{:},1)) Nt]);
A3Data = reshape(A3Data,[size(AData(:,:,id{:},1)) Nt]);

%% Convert to TVMAT
A1 = tvmat(A1Data,T,AIM);
A2 = tvmat(A2Data,T,AIM);
A3 = tvmat(A3Data,T,AIM);