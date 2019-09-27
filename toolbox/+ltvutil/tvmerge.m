function out = tvmerge(varargin)
%% TVMERGE merges the TVMATs on specified horizons

% Total number of TVMATs to be merged
N = nargin;
out = varargin{1};

for i = 2:N
    out = LOCALStitch(out,varargin{i});
end
end

function C = LOCALStitch(A,B)
% Make sure both inputs are finite horizon
ltvutil.verifyFH(A);
ltvutil.verifyFH(B);

% Make sure that Time, InterpolationMethod, SampleTime are consistent
Tf = A.Time(end);
T0 = B.Time(1);
Ainterp = A.InterpolationMethod;
Binterp = B.InterpolationMethod;
ATs = A.Ts;
BTs = B.Ts;
if ~(isequal(Tf,T0) && isequal(Ainterp,Binterp) && isequal(ATs,BTs))
    error('Data, Time, InterpolationMethod and SampleTime must be consistent');
end

% Construct C
Bid = B.Time>T0;
Cinterp = Ainterp;
C = tvmat([squeeze(A.Data);squeeze(B.Data(:,:,Bid))],[A.Time;B.Time(Bid)],Cinterp);
end