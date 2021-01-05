function [isTimeInvariant,varargout] = tvswitchyard(varargin)
% TVSWITCHYARD Utility function for LTV object interaction.
%
% This function checks all inputs for consistency and lifts objects
% to be TVMAT/TVSS. Specifically, the function:
%  1) Raises empty, double, or logical inputs to be a constant TVMAT/TVSS.  
%  2) Checks that all inputs have a consistent Time & InterpolationMethod. 
%  3) If inputs contain both Constant/Varying objects then the Constant
%     objects are mapped to a Varying object on the correct time vector.
%
% Note: Gridded objects combinations with different time vectors are not
% allowed. Use EVALT to convert to a consistent time vector.

%% Get Time and DataTypes
% Assume each input is either a tvss, tvmat, empty, double, or logical. 
% Empty,  double, or logical are lifted to a Constant tvmat/tvss.
nin = nargin;
varargout           = cell(nin,1);
AllTime             = cell(nin,1);
AllIM               = cell(nin,1);
AllisTimeInvariant  = false(nin,1);
AllTs               = zeros(nin,1);
AllTUnit            = cell(nin,1);
for i=1:nin
    % Lift empty/double/logical to a TVMAT or TVSS. This will error out
    % if varargin{i} is an input of unexpected class.
    vi = varargin{i};
    if isa(vi,'double') || isa(vi,'logical')
        vi = tvmat(vi);
    elseif isa(vi,'ss') || isa(vi,'tf') || isa(vi,'zpk')
        vi = tvss(vi);
    elseif isa(vi,'umat') || isa(vi,'ureal') || isa(vi,'ucomplex') || isa(vi,'ucomplexm')
        vi = tvumat(vi);
    elseif isa(vi,'uss') || isa(vi,'ultidyn') || isa(vi,'udyn')
        vi = tvuss(vi);
    end
     
    varargout{i} = vi;
    AllTime{i} = vi.Time;
    AllIM{i} = vi.InterpolationMethod;
    AllTs(i) = vi.Ts;
    AllTUnit{i} = vi.TimeUnit;
    if vi.isTimeInvariant
        AllisTimeInvariant(i) = true;
    end
end

%% Check Constant/Varying
% isTimeInvariant=true if all inputs are constant
isTimeInvariant = all(AllisTimeInvariant);

%% Check Time Vectors
if ~isTimeInvariant
    DTidx = find(AllisTimeInvariant==isTimeInvariant);
    Time = AllTime{ DTidx(1) };
    for i=2:numel(DTidx)
        if ~isequal(Time,AllTime{DTidx(i)} )
            error('Operations require a consistent Time vector.');
        end
    end
end

%% Check InterpolationMethod
if ~isTimeInvariant
    IM = AllIM{ DTidx(1) };
    for i=2:numel(DTidx)
        if ~isequal(IM,AllIM{DTidx(i)} )
            error('Operations require a consistent InterpolationMethod.');
        end
    end
end

%% Check SampleTime
if ~isTimeInvariant
    Ts = AllTs( DTidx(1));
    for i=2:numel(DTidx)
        if ~isequal(Ts,AllTs(DTidx(i)))
            error('Operations require a consistent SampleTime.');
        end
    end
end

%% Check TimeUnit
if ~isTimeInvariant
    TUnit = AllTUnit{ DTidx(1) };
    for i=2:numel(DTidx)
        if ~isequal(TUnit,AllTUnit{DTidx(i)} )
            error('Operations require a consistent TimeUnit.');
        end
    end
end

%% Create output
if ~isTimeInvariant
    for i=1:nin
        if AllisTimeInvariant(i)
            % Lift Constant to TVMAT/TVSS
            Data = varargout{i}.Data;
            nd = ndims(Data)-2;
            Data2 = repmat(Data,[1 1 ones(1,nd) numel(Time)]);            
            if isa(Data,'double') || isa(Data,'logical')
                if isequal(Ts,0)
                    varargout{i} = tvmat(Data2,Time,IM);
                else
                    varargout{i} = tvmat(Data2,Time,Ts);
                end
            elseif isa(Data,'ss') || isa(Data,'tf') || isa(Data,'zpk')
                varargout{i} = tvss(Data2,Time,IM);
            elseif isa(Data,'umat')
                varargout{i} = tvumat(Data2,Time,IM);
            elseif isa(Data,'uss')
                varargout{i} = tvuss(Data2,Time,IM);
            end
        end
    end
end