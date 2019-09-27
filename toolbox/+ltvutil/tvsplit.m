function varargout = tvsplit(A,aTf)
%% TVSPLITEngine splits the TVOBJECT on specified time intervals
%
% [B,C] = TVSPLIT(A,T1) splits the horizon of the TVOBJECT A and
% returns B as A in [T0 T1]
%         C as A in [T1 Tf]
%
% [B,C,D,E,...] = TVSPLIT(A,[T1,T2,T3,...]) splits the horizon of the
% TVOBJECT A and returns B as A in [T0 T1]
%         C as A in [T1 T2]
%         D as A in [T2 T3]
%         E as A in [T3 Tf]
%
% NOTE: Here, [T0,Tf] = getHorizon(A) are assumed to be boundry points,
% input tvmat must be finite horizon and T1, T2, T3 must be sorted such
% that T0 < T1 < T2 < T3 < ... < Tf.
%
% Example:
% t = 0:0.01:10;
% A = tvmat(sin(t),t);
% [B,C,D,E] = tvsplit(A,[3 5 8]);

% Time horizon array for splitting the TVMAT
ltvutil.verifyFH(A);
[T0,Tf] = getHorizon(A);

% Make sure that T0<=t1 and Tf>=tend
saTf = sort(aTf(:));
t1 = saTf(1);
tend = saTf(end);
if ~(T0<=t1 && Tf>=tend)
    error('The new time horizon must be within the original time horizon.')
end
aTf = [T0; saTf; Tf];
nT = length(aTf);
Tcell = cell(nT-1,1);

% Create Horizon Array
for i = 2:nT
    Tcell{i-1} = [aTf(i-1) aTf(i)];
end

% If user included end points then remove those entries
if diff(Tcell{1}) == 0
    Tcell(1) = [];nT = nT - 1;
end
if diff(Tcell{end}) == 0
    Tcell(end) = [];nT = nT - 1;
end

% For Loop
varargout = cell(nT-1,1);
for j = 2:nT
    varargout{j-1} = LOCALtvsplit(A,Tcell{j-1});
end
end

function B = LOCALtvsplit(A,newH)
% Access new horizon being [t1,t2]
t1 = newH(1);
t2 = newH(2);

% Split the time vector inclusive of terminal points
id = (A.Time >= t1 & A.Time <= t2);
requestedTime = A.Time(id);
if requestedTime(end) < t2
    B = evalt(A,[requestedTime;t2]);
else
    B = evalt(A,requestedTime);
end
end