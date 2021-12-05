function H = tvplotEngine(plotfh,varargin)
%% TVPLOTENGINE
% 
% This function plots the time-varying matrix data as per the function
% handle specified
%
% Function handle can be @plot, @stem, @semilogy, @semilogx, @loglog

%% Input Processing
% Replace each TVMAT input by its (reshaped) data
nin = nargin-1;
IgnoreTime = false;

% Convert to strings to chars as strcmp does not work on cell arrays that
% have strings in them.
varargin = controllib.internal.util.hString2Char(varargin);

% Keep track of all the locations of tvmat
idx = zeros(1,nin);
for i=1:nin
    if isa(varargin{i},'tvmat')
        idx(i) = 1;
    end
    if isa(varargin{i},'double')
        idx(i) = 2;
    end
    if isa(varargin{i},'char')
        idx(i) = 3;
    end
end

% Count total number of TVMATs
totalTVMAT = sum(idx==1);

% Preallocate the argIN array
argIN = cell(length(varargin)+totalTVMAT,1);
TU = cell(totalTVMAT,1);
k = 1;

% Decompose the TVMAT into its data and time vector
for i=1:nin
    vi = varargin{i};
    
    % Shift TVMAT Data so that first dimension corresponds to Time
    if isa(vi,'tvmat')
        if ~any(isinf(vi.Time))
            argIN{k} = vi.Time;
            argIN{k+1} = reshapedata(vi);
            TU{i} = vi.TimeUnit;
            k = k + 2;
        else
            error('ltvtools:tvplot:constants',...
                'Infinite Dimensional or Constant TVMAT objects cannot be plotted.')
        end
    else
        argIN{k} = vi;
        k = k + 1;
    end
end

% TVMAT-Data sequence
TVMAT_DATA_SEQ = strfind(idx,[1 2]);
TVMAT_TVMAT_SEQ = strfind(idx,[1 1]);
DATA_TVMAT_SEQ = strfind(idx,[2 1]);

% Remove Time Entry for the above sequences
if ~isempty(TVMAT_TVMAT_SEQ)
    Tvec = cell(numel(TVMAT_TVMAT_SEQ)*2,1);
    k = 1;IgnoreTime = true;
    for i = 1:numel(TVMAT_TVMAT_SEQ)
        % For each sequence entry remove 2 time entry
        SEQ = TVMAT_TVMAT_SEQ(i);
        for j = 1:2
            % Store time vector in Tvec
            Tvec{k} = argIN{SEQ};
            argIN(SEQ) = [];
            
            % Increment k
            k = k + 1;
            SEQ = SEQ + 1;
        end
    end
elseif ~isempty(TVMAT_DATA_SEQ)
    Tvec = cell(numel(TVMAT_DATA_SEQ),1);
    k = 1;IgnoreTime = true;
    for i = 1:numel(TVMAT_DATA_SEQ)
        % Store time vector in Tvec
        Tvec{k} = argIN{TVMAT_DATA_SEQ(i)};
        argIN(TVMAT_DATA_SEQ(i)) = [];
        k = k + 1;
    end
elseif ~isempty(DATA_TVMAT_SEQ)
    Tvec = cell(numel(DATA_TVMAT_SEQ),1);
    k = 1;IgnoreTime = true;
    for i = 1:numel(DATA_TVMAT_SEQ)
        % Store time vector in Tvec
        Tvec{k} = argIN{DATA_TVMAT_SEQ(i)+1};
        argIN(DATA_TVMAT_SEQ(i)+1) = [];
        k = k + 1;
    end
end

% Send argIN to plot function handle
if nargout>0
    H = plotfh(argIN{:});
else
    plotfh(argIN{:});
end

% If IgnoreTime = false then xlabel as Time in TimeUnits
if ~IgnoreTime
    TUu = unique(TU(~cellfun('isempty',TU)));
    if numel(TUu) == 1
        xlabel(sprintf('Time (%s)',TUu{1}));
    else
        xlabel('Time');
    end
end