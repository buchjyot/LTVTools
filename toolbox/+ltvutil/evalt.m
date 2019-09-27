function varargout = evalt(tvflag,varargin)
%% EVALT Engine File
% Main for loop for input arguments
%
% TVFLAG is set to true if return arguments are requested to be in native
% objects

% Input Processing
Tgrid = varargin{end};
Nt = numel(Tgrid);
if tvflag && (~isa(Tgrid,'double')|| Nt<2 || any(diff(Tgrid)< 0))
    error('Last input argument must be a non-decreasing vector of time grid points with length>=2.');
end
nin = nargin-1;
N = nin-1;

% Loop through input arguments
varargout = cell(1,N);
for i = 1:N
    vin = varargin{i};
    varargout{i} = EVALTRouting(vin,Tgrid,tvflag);
end
end

function out = EVALTRouting(in,T,tvflag)
%% EVALTRouting routs the incoming TVOBJECT to respective helper
%
% in     : Incoming TVOBJECT
% T      : Time Grid
% tvfalg : Whether the output should be in native object or array

% Identify class of incoming TVOBJECT and redirect to respective helper
switch class(in)
    case 'tvmat'
        out = LOCALevalTVMAT(in,T,tvflag);
    case 'tvss'
        out = LOCALevalTVSS(in,T,tvflag);
    case 'tvumat'
        out = LOCALevalTVUMAT(in,T,tvflag);
    case 'tvuss'
        out = LOCALevalTVUSS(in,T,tvflag);
end
end

function P = LOCALevalTVSS(G,T,tvflag)
%% LOCALevalTVSS
% Evalute individual state-space matrices on a time grid
[A,B,C,D] = ssdata(G);
At = LOCALevalTVMAT(A,T,tvflag);
Bt = LOCALevalTVMAT(B,T,tvflag);
Ct = LOCALevalTVMAT(C,T,tvflag);
Dt = LOCALevalTVMAT(D,T,tvflag);

% Package up result
if tvflag
    P = tvss(At,Bt,Ct,Dt);
else
    P = ss(At,Bt,Ct,Dt);
end
end

function B = LOCALevalTVUMAT(A,T,tvflag)
%% LOCALevalTVUMAT
% Handle Constant Case
AData = A.Data;
Nd = ndims(AData);
Nt = numel(T);
if A.isTimeInvariant
    BData = repmat(AData,[ones(1,Nd) Nt]);
    if tvflag
        B = tvumat(BData,T);
    else
        B = BData;
    end
    return
else
    error('This is not supported yet.');
end
end

function P = LOCALevalTVUSS(G,T,tvflag)
%% LOCALevalTVUSS

% Use lftdata to extract nominal and uncertain part
[M,Delta] = lftdata(G);

% Interpolate individual elements
M = ltvutil.evalt(tvflag,M,T);

% Wrap back into lft
% XXX: Time Consuming
lftIC = lft(Delta,M);

% Package up result
if tvflag
    P = lftIC;
else
    P = lftIC.Data;
end
end

function B = LOCALevalTVMAT(A,T,tvflag)
%% LOCALevalTVMAT
% Handle Constant Case
AData = A.Data;
Nd = ndims(AData);
Nt = numel(T);
if A.isTimeInvariant
    BData = repmat(AData,[ones(1,Nd) Nt]);
    if tvflag
        B = tvmat(BData,T);
    else
        B = BData;
    end
    return
end

%% Handle Discrete Time-Varying Case
ATime = A.Time;
ATs = A.Ts;

if ~isequal(ATs,0)
    LIA = ismembertol(ATime,T(:));
    if any(LIA) && sum(LIA) == Nt
        BData = AData(:,:,LIA);
    else
        error('Invalid Time input for Discrete-Time TVMAT evalt.');
    end
    % Output
    if tvflag
        B = tvmat(BData,T,A.SampleTime,A.InterpolationMethod,A.TimeUnit);
    else
        B = BData;
    end
    return
end

%% Handle Continuous Time-Varying Case
ATime = A.Time;
AIM = A.InterpolationMethod;

% Pre-process TVMAT Data
ADiff = diff(ATime);
NATime = numel(ATime);
nad = ndims(AData)-3;
if nad==0
    id =  cell(1,0);
else
    id = repmat({':'},1,nad);
end

% Memory Allocation
WarningFlag = false;
mA = size(AData,1);
nA = size(AData,2);
ZEROS = zeros(mA,nA);
BData = zeros( [numel(ZEROS) Nt] );

% Perform Interpolation
switch AIM
    case 'Linear'
        for i=1:Nt
            if T(i)<=ATime(end) && T(i)>=ATime(1)
                % Within the horizon inclusive of boundary points
                [k,alpha] = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
                if alpha==0
                    Bi = AData(:,:,id{:},k);
                else
                    Bi = (1-alpha)*AData(:,:,id{:},k) + ...
                        alpha*AData(:,:,id{:},k+1);
                end
            else
                % Extrapolate to zeros if outside horizon
                Bi = ZEROS;
                WarningFlag = true;
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
        
    case 'Flat'
        for i=1:Nt
            if T(i)<=ATime(end) && T(i)>=ATime(1)
                % Within the horizon inclusive of boundry points
                k = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
                Bi = AData(:,:,id{:},k);
            else
                % Extrapolate to zeros if outside horizon
                Bi = ZEROS;
                WarningFlag = true;
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
        
    case 'Nearest'
        for i=1:Nt
            if T(i)<=ATime(end) && T(i)>=ATime(1)
                [k,alpha] = LOCALfindslotalpha(NATime,ATime,T(i),ADiff);
                if alpha>=0.5
                    k=k+1;
                end
                % Within the horizon inclusive of boundry points
                Bi = AData(:,:,id{:},k);
            else
                % Extrapolate to zeros if outside horizon
                Bi = ZEROS;
                WarningFlag = true;
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(Bi) Nt]);
        
    case 'Spline'
        % Currently SplineData is empty and it computed in every call
        % to EVALT. A syntax to precompute Spline Data is:
        %    A = getSplineData(A);
        % This allows for more efficient, repeated calls to EVALT.
        ASplineData = A.SplineData;
        if isempty(ASplineData)
            ASplineData = getSplineData(AData,ATime);
        end
        
        for i=1:Nt
            if T(i)<ATime(1) || T(i)>ATime(end)
                % Extrapolate to zeros if outside horizon
                Bi = ZEROS;
                WarningFlag = true;
            else
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
                    Bi = ASplineData(aidx(1),:);
                else
                    Bi = [1 dt dt^2 dt^3]*ASplineData(aidx,:);
                end
            end
            BData(:,i) = Bi(:);
        end
        BData = reshape(BData,[size(AData(:,:,id{:},1)) Nt]);
end

%% Convert output to a TVMAT if tvflag specified
if tvflag
    B = tvmat(BData,T,AIM);
else
    B = BData;
end

% Throw warning if time grid points are outside horizon
if WarningFlag
    warning('ltvtools:evalt:extrapolate',...
        'Requested time grid points are outside the horizon, extrapolating to zeros.');
end
end

function [k,alpha] = LOCALfindslotalpha(N,vec,val,dvec)
%% LOCALfindslotalpha
% N integer
% vec 1-by-N (or N-by-1), sorted
% val 1-by-1
% dvec = diff(vec)

k = max(find(val>=vec));   %#ok<MXFND> % don't follow advice - it is slower.
if ~isempty(k)
    if k<N
        alpha = (val - vec(k))/dvec(k);
    else
        alpha = 0;
    end
else
    k = 1;
    alpha = 0;
end
end