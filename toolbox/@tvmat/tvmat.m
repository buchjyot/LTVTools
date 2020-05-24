%% TVMAT Create a time-varying matrix
%
%   M = tvmat(Data,Time) creates a time-varying matrix. Time is an Nt-by-1
%   vector of increasing time values.  If Data is 1-by-Nt then M is a
%   time-varying scalar with value Data(i) at Time(i). If Data is N1-by-Nt
%   then M is an N1-by-1 time-varying vector with value Data(:,i) at
%   Time(i). In general, if Data is an N1-by-N2-by-...-by-Nad-by-Nt array
%   then M is an N1-by-N2-by-...-by-Nad time-varying array with value
%   Data(:,...,:,i) at Time(i).
%
%   M = tvmat(Data,Time,InterpolationMethod) specifies the method to
%   interpolate the data between the specified time grid points.  The
%   method can be 'Linear', 'Flat', 'Nearest', or 'Spline'.
%   Default is InterpolationMethod = 'Linear'.
%
%   M = tvmat(Data) generates a constant matrix.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 matrix defined on a 1-dimensional grid
%
%   M1 = tvmat(5);
%   Time = linspace(0,5,10)';
%   Data = randn(1,1,10);
%   M2 = tvmat(Data,Time);
%   M3 = M1+M2;
%   plot(M2.Time,M2,'b',M3.Time,M3,'r--');
%
% See also:

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd}) tvmat < matlab.mixin.CustomDisplay
    % Class definition for time-varying matrix
    
    properties (SetAccess = public)
        Data = [];
        Time = [-inf,inf];
        InterpolationMethod = 'Linear';
        % XXX: Add flag to compute spline data?  (We don't want to compute
        % this every time the constructor is called as this causes all
        % operations to be slow.  However, we also need some method to say
        % that we want to compute it otherwise it will be computed
        % everytime EVALT is called)
        % Store offset data? (e.g. as returned by linearize)
    end
    
    properties (SetAccess = immutable)
        % Can not be changed by user
        isTimeInvariant = false;
    end
    
    properties (SetAccess = private)
        % Note: This class has a subsasgn method which, in general, by-
        % passes "SetAccess = private". The private SetAccess is actually
        % enforced by logic in subsasgn.
        SplineData = [];
    end
    
    properties
        % Following Properties are not directly displayed but are
        % customized, see protected methods of this class or use "get"
        % command to see them
        Ts = 0;
        TimeUnit = 'seconds';
        Name = '';
    end
    
    methods
        %% Constructor
        % XXX Add undocumented flag to skip error checking?
        function obj = tvmat(Data,Time,InterpolationMethod,varargin)
            nin = nargin;
            switch nin
                case 0
                    % Return default object
                    return
                case 1
                    if isa(Data,'tvmat')
                        obj = Data;
                        return
                    else
                        obj.Data = Data;
                        obj.isTimeInvariant = true;
                    end
                case 2
                    % By Default Continuous-Time
                    obj.Data = Data;
                    obj.Time = Time(:);
                case 3
                    obj.Data = Data;
                    obj.Time = Time(:);
                    obj.InterpolationMethod = InterpolationMethod;
                case 4
                    obj.Data = Data;
                    obj.Time = Time(:);
                    obj.InterpolationMethod = InterpolationMethod;
                    obj.Ts = varargin{1};
                case 5
                    obj.Data = Data;
                    obj.Time = Time(:);
                    obj.InterpolationMethod = InterpolationMethod;
                    obj.Ts = varargin{1};
                    obj.TimeUnit = varargin{2};
            end
            
            % User want us to compute the SampleTime
            if isempty(obj.Ts) || isempty(obj.InterpolationMethod)
                diffTime = diff(Time(:));
                obj.Ts = uniquetol(diffTime,sqrt(eps));
            end
            
            Nt = numel(obj.Time);
            if ~obj.isTimeInvariant
                % Code below allows Data to be Nt-by-Nd or Nd-by-Nt.
                % It reshapes Data to Nd-by-1-by-Nt.
                if ismatrix(Data)
                    szData = size(Data);
                    if szData(2)==Nt
                        obj.Data = reshape(Data,[szData(1) 1 Nt]);
                    elseif szData(1)==Nt
                        obj.Data = reshape(Data',[szData(2) 1 Nt]);
                    end
                end
            end
            
            % Squeeze extra singleton dimensions
            if ~obj.isTimeInvariant && ~all(size(Data)==0)
                szData = size(obj.Data);
                szData(end) = 1;
                idx = find( diff(szData)~= 0 );
                Nd = max(2, numel(idx) );
                obj.Data = reshape(obj.Data, [szData(1:Nd) Nt]);
            end
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        %% ISVALID
        % isvalid: Determine if TVMAT object has valid data
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if TVMAT object is valid.
            errstr = [];
            pflag = 1;
            
            % Check Data Class
            objData = obj.Data;
            if ~( isa(objData,'double') || isa(objData,'logical')  )
                pflag = 0;
                errstr = 'Data must be a double or logical.';
                return
            end
            
            % Check Time
            objTime = obj.Time;
            Nt = numel(objTime);
            if ~isa(objTime,'double') || Nt<2 || any( diff(objTime)< 0 )
                pflag = 0;
                errstr = ['Time vector must be non-decreasing vector ' ...
                    'of doubles with length>=2.'];
                return
            end
            
            % Check SampleTime
            objTs = obj.Ts;
            if ~(isnumeric(objTs) && isscalar(objTs) && isreal(objTs) && isfinite(objTs))
                error('The value of the "SampleTime" property must be a real number.');
            elseif objTs<0 && objTs~=-1
                error('The value of the "SampleTime" property must be 0, a positive scalar, or -1 to mean unspecified for constant TVMATs.');
            end
            if ~isequal(objTs,0)
                error('Discrete-Time TVMAT are not supported yet.');
            end
            
            % Check Time/Data Consistency
            szData = size(obj.Data);
            if any(isnan(objTime))
                pflag = 0;
                errstr = 'Time vector must not contain any NaNs.';
            end
            if obj.isTimeInvariant && ~all(objTime==[-inf inf])
                pflag = 0;
                errstr = 'Time vector must be [-inf inf] for constant TVMATs.';
            end
            if ~obj.isTimeInvariant && ~all(size(objData)==0) && ...
                    ~(numel(szData)>2 && szData(end)==Nt )
                pflag = 0;
                errstr = 'Dimensions of Time and Data are incompatible.';
            end
            if ~isequal(objTs,0) && ~isequal(objTs,-1)
                % Discrete-Time TVMAT
                diffTime = diff(objTime);
                
                if ~isequal(numel(uniquetol(diffTime,sqrt(eps))),1)
                    pflag = 0;
                    errstr = 'Discrete-Time TVMAT must be regularly sampled.';
                elseif ~isequal(numel(uniquetol([uniquetol(diffTime,sqrt(eps)),...
                        obj.Ts],sqrt(eps))),1)
                    pflag = 0;
                    errstr = 'SampleTime Property must be consistent with Time.';
                end
            end
            
            % Check IM
            IM = obj.InterpolationMethod;
            if ~isequal(objTs,0)
                % Discrete-Time TVMAT IM should be empty
                if ~isempty(IM)
                    pflag = 0;
                    errstr = 'InterpolationMethod must be empty for Discrete-Time TVMAT.';
                end
            elseif ~any(strcmpi(IM,{'Flat'; 'Nearest';'Linear'; 'Spline'}) )
                pflag = 0;
                errstr = ['InterpolationMethod must be "Flat", '...
                    '"Nearest", "Linear", or "Spline" for Continuous-Time TVMAT'];
            end
            
            if pflag==0 && nargout==0
                error(errstr);
            end
        end
        
        %% GET Method
        function V = get(obj,Property)
            % GET  Access/query property values.
            %
            %   VALUE = GET(OBJ,'PropertyName') returns the value of the
            %   specified property of the input/output model M. An
            %   equivalent syntax is VALUE = M.PropertyName.
            %
            %   VALUES = GET(OBJ,{'Prop1','Prop2',...}) returns the values
            %   of several properties at once in a cell array VALUES.
            %
            %   GET(OBJ) displays all properties of M and their values.
            %
            %   S = GET(OBJ) returns a structure whose field names and
            %   values are the property names and values of M.
            
            narginchk(1,2);
            if nargin==2
                % GET(M,'Property') or GET(M,{'Prop1','Prop2',...})
                CharProp = ischar(Property);
                if CharProp
                    Property = {Property};
                else
                    try
                        Property = cellstr(Property);
                    catch ME
                        error(message('Control:ltiobject:get1'))
                    end
                end
                AllPublicProps = ltipack.allprops(obj);
                ClassName = class(obj);
                
                % Get all public properties
                Nq = numel(Property);
                V = cell(1,Nq);
                for i=1:Nq
                    try
                        V{i} = obj.(ltipack.matchProperty(Property{i},AllPublicProps,ClassName));
                    catch E
                        throw(E)
                    end
                end
                
                % Strip cell header if PROPERTY was a string
                if CharProp
                    V = V{1};
                end
            else
                % Construct struct of field values
                propList = [];
                allProp = properties(obj);
                for i = 1:length(allProp)
                    propList.(allProp{i}) = obj.(allProp{i});
                end
                if nargout
                    V = propList;
                else
                    disp(propList)
                end
            end
        end
        
        %% SET Methods
        function obj = set.Ts(obj,V)
            obj.Ts = V;
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        function obj = set.TimeUnit(obj,V)
            % SET function for TimeUnit property
            if strcmp(V,'')
                TU = 'seconds';  % for SITB backward compatibility
            else
                TU = ltipack.matchKey(V,ltipack.getValidTimeUnits());
                if isempty(TU)
                    error(message('Control:ltiobject:setTimeUnit'))
                end
            end
            if ~strcmp(TU,obj.TimeUnit)
                obj.TimeUnit = TU;
            end
        end
        
        function obj = set.Name(obj,V)
            if ischar(V)
                obj.Name = V;
            else
                error('Property "Name" must be set to a character vector, for example, ''abc''')
            end
        end
        
        function out = set(obj,varargin)
            % SET  Modifies values of TVMAT properties.
            %
            %   OBJ = SET(OBJ,'Property',VALUE) sets the property with name 'Property'
            %   to the value VALUE. This is equivalent to SYS.Property = VALUE.
            %   SYS can be any dynamic system object.
            %
            %   OBJ = SET(OBJ,'Property1',Value1,'Property2',Value2,...) sets multiple
            %   property values in a single command.
            %
            %   OBJ = SET(OBJ,'Property1',Value1,...) returns the modified
            %   system SYSOUT.
            %
            %   SET(OBJ,'Property') displays information about valid values for
            %   the specified property of SYS.
            %
            %   See also GET.
            
            ni = nargin;
            no = nargout;
            ClassName = class(obj);
            
            if ni<3
                % Informative syntax
                if no>0
                    error(message('Control:general:OptionHelper1'))
                else
                    error('The "set" command requires at least three input arguments.');
                end
            else
                % SET(SYS,'Prop1',Value1, ...)
                % Set property values, deferring the consistency checks until all
                % properties have been set
                ltipack.mustBeNameValuePairs(varargin)
                PublicProps = ltipack.allprops(obj);
                for ct=1:2:ni-1
                    obj.(ltipack.matchProperty(varargin{ct},PublicProps,ClassName)) = varargin{ct+1};
                end
                
                if no>0
                    out = obj;
                else
                    % Use ASSIGNIN to update in place
                    sysname = inputname(1);
                    if isempty(sysname)
                        error(message('Control:ltiobject:setLTI5'))
                    end
                    assignin('caller',sysname,obj)
                end
            end
        end
        
        %% Basic Functions
        
        % XXX TODO
        % numArgumentsFromSubscript
        
        function out = double(obj)
            out = obj.Data;
        end
        
        function out = isct(obj)
            % Returns true if SampleTime is 0
            out = isequal(obj.Ts,0);
        end
        
        function out = isdt(obj)
            % Returns true is SampleTime is nonzero
            out = ~isequal(obj.Ts,0);
        end
        
        function out = isfh(obj)
            % Returns true if object is finite horizon
            out = ~any(isinf(obj.Time));
        end
        
        function out = isUncertain(obj) %#ok<MANU>
            % Returns true if Uncertain object
            out = false;
        end
        
        function [T0,Tf] = getHorizon(obj)
            % Returns the horizon of the object
            T0 = obj.Time(1);
            Tf = obj.Time(end);
        end
        
        function out = end(obj,slot,nslots)
            sz = size(obj);
            if slot < nslots
                out = sz(slot);
            else
                out = prod(sz(slot:end));
            end
        end
        
        function out = length(obj)
            out = size(obj);
            out = max(out);
        end
        
        function [Y,I] = min(X,varargin)
            if nargin>=2 %&& isa(varargin{1},'tvmat')
                Y = varargin{1};
                [Y,I] = tvbinop(@min,X,Y,varargin{2:end});
            else
                [Y,I] = tvunop(@min,X,varargin{:});
            end
        end
        
        function [Y,I] = max(X,varargin)
            if nargin>=2 %&& isa(varargin{1},'tvmat')
                Y = varargin{1};
                [Y,I] = tvbinop(@max,X,Y,varargin{2:end});
            else
                [Y,I] = tvunop(@max,X,varargin{:});
            end
        end
        
        function out = tvmax(obj)
            out = max(obj.Data(:));
        end
        
        function out = tvmin(obj)
            out = min(obj.Data(:));
        end
        
        function out = isempty(obj)
            %szo = size(obj); out = any(~szo);
            out = isempty(obj.Data);
        end
        
        function out = iscolumn(obj)
            out = size(obj,2)==1;
        end
        
        function out = isrow(obj)
            out = size(obj,1)==1;
        end
        
        function out = isscalar(obj)
            szo = size(obj);
            out = (szo(1)==1 && szo(2)==1);
        end
        
        function out = isvector(obj)
            szo = size(obj);
            out = (szo(1)==1 || szo(2)==1);
        end
        
        function out = ndims(obj)
            if obj.isTimeInvariant
                out = ndims( obj.Data );
            else
                out = ndims( obj.Data ) - 1;
            end
        end
        
        function varargout = size(obj,varargin)
            if obj.isTimeInvariant
                [varargout{1:max(nargout,1)}] = size( obj.Data, varargin{:});
            else
                objData = obj.Data;
                nd = ndims(objData)-3;
                id = repmat({':'},1,nd);
                [varargout{1:max(nargout,1)}] = size( objData(:,:,id{:},1), varargin{:});
            end
        end
        
        function varargout = tv2fh(varargin)
            % Returns function handle for a corresponding time-varying object
            % FH = tv2fh(TVMAT)
            nin = nargin;
            varargout = cell(nin,1);
            for i = 1:nin
                varargout{i} = ltvutil.tv2fh(varargin{i});
            end
        end
        
        function varargout = tv2ts(varargin)
            % Interface to timeseries objects
            % TS = tv2ts(TVMAT)
            nin = nargin;
            varargout = cell(nin,1);
            for i = 1:nin
                M = varargin{i};
                MATName = M.Name;
                ts = timeseries(M.Data,M.Time,'Name',MATName);
                ts.TimeInfo.Units = M.TimeUnit;
                if contains(M.InterpolationMethod,{'Linear','ZOH'})
                    varargout{i} = setinterpmethod(ts,M.InterpolationMethod);
                else
                    warning('ltvtools:tvmat:convert2ts',...
                        'Timeseries objects only support Linear or ZOH interpolation method.');
                    varargout{i} = ts;
                end
            end
        end
        
        function varargout = evalt(varargin)
            % EVALT Evaluate a TVMAT(s) at specified time(s)
            %
            % B = EVALT(A,T) evaluates a TVMAT A on the vector of times T
            % using the interpolation method specified by
            % A.InterpolationMethod.  B is returned as a TVMAT.
            %
            % [C,D,...] = EVALT(A,B,...,T) evaluates a TVMAT A and TVMAT B
            % etc. on the vector of times T. [C,D,...] is returned as
            % TVMATs.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(true,varargin{:});
        end
        
        function varargout = tvsubs(varargin)
            % TVSUBS Evaluate a TVMAT(s) at specified time(s)
            %
            % B = TVSUBS(A,T) evaluates a TVMAT A on the vector of times T
            % using the interpolation method specified by
            % A.InterpolationMethod.  B is returned as a double.
            %
            % [C,D,...] = TVSUBS(A,B,...,T) evaluates a TVMAT A and TVMAT B
            % etc. on the vector of times T. [C,D,...] is returned as
            % double.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(false,varargin{:});
        end
        
        function varargout = tvsplit(A,aTf)
            % TVSPLIT splits the TVMAT on specified time intervals
            %
            % [B,C] = TVSPLIT(A,T1) splits the horizon of the tvmat A and
            % returns B as A in [T0 T1]
            %         C as A in [T1 Tf]
            %
            % [B,C,D,E,...] = TVSPLIT(A,[T1,T2,T3,...]) splits the horizon
            % of the tvmat A and returns B as A in [T0 T1]
            %         C as A in [T1 T2]
            %         D as A in [T2 T3]
            %         E as A in [T3 Tf]
            %
            % NOTE: Here, [T0,Tf] = getHorizon(A) are assumed to be boundry
            % points, input tvmat must be finite horizon and T1, T2, T3
            % must be sorted such that T0 < T1 < T2 < T3 < ... < Tf.
            %
            % Example:
            % t = 0:0.01:10;
            % A = tvmat(sin(t),t);
            % [B,C,D,E] = tvsplit(A,[3 5 8]);
            [varargout{1:max(nargout,1)}] = ltvutil.tvsplit(A,aTf);
        end
        
        function varargout = tvmerge(varargin)
            % TVMERGE merges the TVMATs on specified horizons
            [varargout{1:max(nargout,1)}] = ltvutil.tvmerge(varargin{:});
        end
        
        %% Unary Operations: Element-by-Element
        function out = abs(obj)
            out = tvunopebe(@abs,obj);
        end
        function out = acos(obj)
            out = tvunopebe(@acos,obj);
        end
        function out = acosd(obj)
            out = tvunopebe(@acosd,obj);
        end
        function out = acosh(obj)
            out = tvunopebe(@acosh,obj);
        end
        function out = acot(obj)
            out = tvunopebe(@acot,obj);
        end
        function out = acotd(obj)
            out = tvunopebe(@acotd,obj);
        end
        function out = acoth(obj)
            out = tvunopebe(@acoth,obj);
        end
        function out = acsc(obj)
            out = tvunopebe(@acsc,obj);
        end
        function out = acscd(obj)
            out = tvunopebe(@acscd,obj);
        end
        function out = acsch(obj)
            out = tvunopebe(@acsch,obj);
        end
        function out = angle(obj)
            out = tvunopebe(@angle,obj);
        end
        function out = asec(obj)
            out = tvunopebe(@asec,obj);
        end
        function out = asecd(obj)
            out = tvunopebe(@asecd,obj);
        end
        function out = asech(obj)
            out = tvunopebe(@asech,obj);
        end
        function out = asin(obj)
            out = tvunopebe(@asin,obj);
        end
        function out = asind(obj)
            out = tvunopebe(@asind,obj);
        end
        function out = asinh(obj)
            out = tvunopebe(@asinh,obj);
        end
        function out = atan(obj)
            out = tvunopebe(@atan,obj);
        end
        function out = atand(obj)
            out = tvunopebe(@atand,obj);
        end
        function out = atanh(obj)
            out = tvunopebe(@atanh,obj);
        end
        function out = ceil(obj)
            out = tvunopebe(@ceil,obj);
        end
        function out = conj(obj)
            out = tvunopebe(@conj,obj);
        end
        function out = cos(obj)
            out = tvunopebe(@cos,obj);
        end
        function out = cosd(obj)
            out = tvunopebe(@cosd,obj);
        end
        function out = cosh(obj)
            out = tvunopebe(@cosh,obj);
        end
        function out = cot(obj)
            out = tvunopebe(@cot,obj);
        end
        function out = cotd(obj)
            out = tvunopebe(@cotd,obj);
        end
        function out = coth(obj)
            out = tvunopebe(@coth,obj);
        end
        function out = csc(obj)
            out = tvunopebe(@csc,obj);
        end
        function out = cscd(obj)
            out = tvunopebe(@cscd,obj);
        end
        function out = csch(obj)
            out = tvunopebe(@csch,obj);
        end
        function out = exp(obj)
            out = tvunopebe(@exp,obj);
        end
        function out = fix(obj)
            out = tvunopebe(@fix,obj);
        end
        function out = fliplr(obj)
            out = tvunopebe(@flip,2,obj);
        end
        function out = flipud(obj)
            out = tvunopebe(@flip,1,obj);
        end
        function out = floor(obj)
            out = tvunopebe(@floor,obj);
        end
        function out = imag(obj)
            out = tvunopebe(@imag,obj);
        end
        function out = isfinite(obj)
            out = tvunopebe(@isfinite,obj);
        end
        function out = isinf(obj)
            out = tvunopebe(@isinf,obj);
        end
        function out = isnan(obj)
            out = tvunopebe(@isnan,obj);
        end
        function out = log(obj)
            out = tvunopebe(@log,obj);
        end
        function out = log10(obj)
            out = tvunopebe(@log10,obj);
        end
        function out = log2(obj)
            out = tvunopebe(@log2,obj);
        end
        function out = not(obj)
            out = tvunopebe(@not,obj);
        end
        function out = real(obj)
            out = tvunopebe(@real,obj);
        end
        function out = round(obj)
            out = tvunopebe(@round,obj);
        end
        function out = sec(obj)
            out = tvunopebe(@sec,obj);
        end
        function out = secd(obj)
            out = tvunopebe(@secd,obj);
        end
        function out = sech(obj)
            out = tvunopebe(@sech,obj);
        end
        function out = sign(obj)
            out = tvunopebe(@sign,obj);
        end
        function out = sin(obj)
            out = tvunopebe(@sin,obj);
        end
        function out = sind(obj)
            out = tvunopebe(@sind,obj);
        end
        function out = sinh(obj)
            out = tvunopebe(@sinh,obj);
        end
        function out = sqrt(obj)
            out = tvunopebe(@sqrt,obj);
        end
        function out = tan(obj)
            out = tvunopebe(@tan,obj);
        end
        function out = tand(obj)
            out = tvunopebe(@tand,obj);
        end
        function out = tanh(obj)
            out = tvunopebe(@tanh,obj);
        end
        function out = uminus(obj)
            out = tvunopebe(@uminus,obj);
        end
        function out = unwrap(varargin)
            out = tvunopebe(@unwrap,varargin{:});
        end
        function out = uplus(obj)
            out = tvunopebe(@uplus,obj);
        end
        function out = rad2deg(obj)
            out = tvunopebe(@rad2deg,obj);
        end
        function out = deg2rad(obj)
            out = tvunopebe(@deg2rad,obj);
        end
        
        %% Unary Matrix Operations: For-loop over time
        function out = cond(obj,varargin)
            out = tvunop(@cond,obj,varargin{:});
        end
        function out = ctranspose(obj)
            out = tvunop(@ctranspose,obj);
        end
        function out = det(obj)
            out = tvunop(@det,obj);
        end
        function out = diag(obj,varargin)
            out = tvunop(@diag,obj,varargin{:});
        end
        function out = expm(obj)
            out = tvunop(@expm,obj);
        end
        function out = inv(obj)
            out = tvunop(@inv,obj);
        end
        function out = mean(obj,varargin)
            out = tvunop(@mean,obj,varargin{:});
        end
        function out = median(obj,varargin)
            out = tvunop(@median,obj,varargin{:});
        end
        function out = norm(obj,varargin)
            out = tvunop(@norm,obj,varargin{:});
        end
        function out = null(obj,varargin)
            out = tvunop(@null,obj,varargin{:});
        end
        function out = permute(obj,varargin)
            out = tvunop(@permute,obj,varargin{:});
        end
        function out = pinv(obj,varargin)
            out = tvunop(@pinv,obj,varargin{:});
        end
        function out = repmat(obj,varargin)
            out = tvunop(@repmat,obj,varargin{:});
        end
        function out = rank(obj,varargin)
            out = tvunop(@rank,obj,varargin{:});
        end
        function out = reshape(obj,varargin)
            out = tvunop(@reshape,obj,varargin{:});
        end
        function [out,t] = reshapedata(obj)
            % Returns reshaped data matrix and a corresponding time vector
            
            % Check Finite Horizon
            ltvutil.verifyFH(obj);
            objData = obj.Data;
            nD = ndims(objData);
            
            % Allocate Memory
            t = obj.Time;
            [m,n,Nt] = size(objData);
            out = zeros(Nt,m*n);
            
            % Main Switch-Case Statements
            switch nD
                case 1
                    out = reshape(objData,size(out));
                case 2
                    out = reshape(objData,Nt,prod(m*n));
                case 3
                    out = reshape(permute(objData,[2 1 3]),prod(m*n),[])';
                otherwise
                    error('TVDATA is not supported for any higher than 3-D Data array.');
            end
        end
        function out = rcond(obj)
            out = tvunop(@rcond,obj);
        end
        function out = squeeze(obj)
            out = tvunop(@squeeze,obj);
        end
        function out = trace(obj)
            out = tvunop(@trace,obj);
        end
        function out = transpose(obj)
            out = tvunop(@transpose,obj);
        end
        
        %% Unary Matrix Operations: For-loop + Variable Outputs
        function varargout = chol(obj,varargin)
            [varargout{1:max(nargout,1)}] = tvunop(@chol,obj,varargin{:});
        end
        function varargout = ldl(obj,varargin)
            [varargout{1:max(nargout,1)}] = tvunop(@ldl,obj,varargin{:});
        end
        function varargout = damp(obj)
            [varargout{1:max(nargout,1)}] = tvunop(@damp,obj);
        end
        function varargout = logm(obj)
            [varargout{1:max(nargout,1)}] = tvunop(@logm,obj);
        end
        function varargout = schur(obj,varargin)
            [varargout{1:max(nargout,1)}] = tvunop(@schur,obj,varargin{:});
        end
        function varargout = sqrtm(obj)
            [varargout{1:max(nargout,1)}] = tvunop(@sqrtm,obj);
        end
        function varargout = sort(obj,varargin)
            [varargout{1:max(nargout,1)}] = tvunop(@sort,obj,varargin{:});
        end
        
        %% Binary Matrix Operations : Single Output
        function out = and(A,B)
            out = tvbinop(@and,A,B);
        end
        function out = atan2(A,B)
            out = tvbinop(@atan2,A,B);
        end
        function out = atan2d(A,B)
            out = tvbinop(@atan2d,A,B);
        end
        function out = eq(A,B)
            out = tvbinop(@eq,A,B);
        end
        function out = ge(A,B)
            out = tvbinop(@ge,A,B);
        end
        function out = gt(A,B)
            out = tvbinop(@gt,A,B);
        end
        function out = ldivide(A,B)
            out = tvbinop(@ldivide,A,B);
        end
        function out = le(A,B)
            out = tvbinop(@le,A,B);
        end
        function out = lft(A,B,varargin)
            out = tvbinop(@lft,A,B,varargin{:});
        end
        function out = lt(A,B)
            out = tvbinop(@lt,A,B);
        end
        function out = minus(A,B)
            out = tvbinop(@minus,A,B);
        end
        function out = mldivide(A,B)
            out = tvbinop(@mldivide,A,B);
        end
        function out = mpower(A,B)
            out = tvbinop(@mpower,A,B);
        end
        function out = mrdivide(A,B)
            out = tvbinop(@mrdivide,A,B);
        end
        function out = mtimes(A,B)
            out = tvbinop(@mtimes,A,B);
        end
        function out = ne(A,B)
            out = tvbinop(@ne,A,B);
        end
        function out = or(A,B)
            out = tvbinop(@or,A,B);
        end
        function out = plus(A,B)
            out = tvbinop(@plus,A,B);
        end
        function out = power(A,B)
            out = tvbinop(@power,A,B);
        end
        function out = rdivide(A,B)
            out = tvbinop(@rdivide,A,B);
        end
        function out = times(A,B)
            out = tvbinop(@times,A,B);
        end
        function out = xor(A,B)
            out = tvbinop(@xor,A,B);
        end
        
        %% N-ary Matrix Operations
        function varargout = care(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@care,varargin{:});
        end
        function varargout = dare(varargin)
            [varargout{1:max(nargout,1)}]= tvnop(@dare,varargin{:});
        end
        function varargout = dlqr(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@dlqr,varargin{:});
        end
        function  out = dlyap(varargin)
            out = tvnop(@dlyap,varargin{:});
        end
        function  out = dlyapchol(varargin)
            out = tvnop(@dlyapchol,varargin{:});
        end
        function varargout = eig(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@eig,varargin{:});
        end
        function varargout = lqr(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@lqr,varargin{:});
        end
        function varargout = lu(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@lu,varargin{:});
        end
        function  out = lyap(varargin)
            out = tvnop(@lyap,varargin{:});
        end
        function  out = lyapchol(varargin)
            out = tvnop(@lyapchol,varargin{:});
        end
        function varargout = qr(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@qr,varargin{:});
        end
        function varargout = svd(varargin)
            [varargout{1:max(nargout,1)}] = tvnop(@svd,varargin{:});
        end
        function  out = sylvester(varargin)
            out = tvnop(@sylvester,varargin{:});
        end
        function out = tril(varargin)
            out = tvnop(@tril,varargin{:});
        end
        function out = triu(varargin)
            out = tvnop(@triu,varargin{:});
        end
        
        %% Recursive Binary Operations
        function out = blkdiag(varargin)
            out = tvrecop(@blkdiag,varargin{:});
        end
        function out = horzcat(varargin)
            out = tvrecop(@horzcat,varargin{:});
        end
        function out = vertcat(varargin)
            out = tvrecop(@vertcat,varargin{:});
        end
        
    end % end of methods
    
    %% DISPLAY Methods
    methods (Access = protected)
        
        %% GetHeader
        function header = getHeader(obj)
            % Identify Type of TVMAT
            if isequal(obj.Ts,0)
                header = ' Continuous-time time-varying matrix';
            elseif ~isequal(obj.Ts,0) && ~isequal(obj.Ts,-1) && ~isempty(obj)
                header = ' Discrete-time time-varying matrix';
            end
            
            % Merge with size
            [Nr,Nc] = size(obj);
            header = [header sprintf(' with %d rows and %d columns.',Nr,Nc)];
            
            % Format
            format = get(0,'FormatSpacing');
            switch format
                case 'compact'
                    header = [newline header newline newline];
                case 'loose'
                    header = [header newline];
            end
        end
        
        %% GetFooter
        function footer = getFooter(obj)             %#ok<MANU>
            footer = ' Use "get" method to see all the properties.';
            format = get(0,'FormatSpacing');
            switch format
                case 'compact'
                    footer = [newline footer newline newline];
                case 'loose'
                    footer = [footer newline newline];
            end
        end
        
        %% Display Empty Object
        function displayEmptyObject(obj)
            displayGeneral(obj);
        end
        
        %% Display Scalar Object
        function displayScalarObject(obj)
            displayGeneral(obj);
        end
        
        %% Display NonScalar Object
        function displayNonScalarObject(obj)
            displayGeneral(obj);
        end
        
        %% Display General
        function displayGeneral(obj)
            if isempty(obj)
                fprintf(newline);
                fprintf(' %s','Empty time-varying matrix.');
                fprintf(newline);fprintf(newline);
            else
                % Header
                fprintf('%s',getHeader(obj));
                
                % Print relevant details
                propList = [];
                if ~isempty(obj.Name)
                    propList.Name = obj.Name;
                end
                propList.Data = obj.Data;
                propList.Horizon = [obj.Time(1),obj.Time(end)];
                if obj.Ts~=0
                    propList.SampleTime = obj.Ts;
                end
                propList.InterpolationMethod = obj.InterpolationMethod;
                propgrp = matlab.mixin.util.PropertyGroup(propList);
                matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgrp);
                
                % Footer
                fprintf('%s',getFooter(obj));
            end
        end
    end
end % end of classdef