%% TVUMAT Create an uncertain time-varying matrix (gridded umat)
%
%   M = tvumat(Data,Time) creates an uncertain time-varying matrix. Time is
%   an Nt-by-1 vector of increasing time values.  If Data is 1-by-Nt then M
%   is a time-varying scalar with value Data(i) at Time(i). If Data is
%   N1-by-Nt then M is an N1-by-1 time-varying vector with value Data(:,i)
%   at Time(i). In general, if Data is an N1-by-N2-by-...-by-Nad-by-Nt
%   array then M is an N1-by-N2-by-...-by-Nad time-varying array with value
%   Data(:,...,:,i) at Time(i).
%
%   M = tvumat(Data,Time,InterpolationMethod) specifies the method to
%   interpolate the data between the specified time grid points.  The
%   method can be 'Linear', 'Flat', 'Nearest', or 'Spline'. Default is
%   InterpolationMethod = 'Linear'.
%
%   M = tvumat(Data) generates a constant uncertain matrix.
%
%   EXAMPLE: (CUT/PASTE)
%   % Create a 2-by-2 uncertain time-varying matrix
%   % defined on a 1-dimensional grid
%
%   p = ureal('p',1);
%   Data = [1 p; p 1];
%   Time = linspace(0,5,10);
%   uM = evalt(tvumat(Data),Time);
%
% See also:

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd,?tvmat,?tvss}) tvumat < matlab.mixin.CustomDisplay
    % Class definition for time-varying umat
    
    properties
        Data = umat;
        Time = [-inf,inf];
        InterpolationMethod = 'Linear';
    end
    
    properties (SetAccess = immutable)
        isTimeInvariant = false;
    end
    
    properties (SetAccess = private)
        % Note: This class has a subsasgn method which, in general, by-
        % passes "SetAccess = private". The private SetAccess is actually
        % enforced by logic in subsasgn.
        SplineData = [];
    end
    
    properties (Dependent)
        % Following properties are dependent on underlying umat array object
        % High-Level Summary of Low-Level umat array
        NominalValue
        Uncertainty
        SamplingGrid
        Name
    end
    
    properties
        % Following Properties are not directly displayed but are
        % customized, see protected methods of this class or use "get"
        % command to see them
        Ts = 0;
        TimeUnit = 'seconds';
    end
    
    methods
        %% Constructor
        % XXX Add undocumented flag to skip error checking?
        function obj = tvumat(Data,Time,InterpolationMethod,varargin)
            nin = nargin;
            switch nin
                case 0
                    % Return default object
                    return
                case 1
                    if isa(Data,'tvumat')
                        obj = Data;
                        return
                    else
                        obj.Data = Data;
                        obj.isTimeInvariant = true;
                    end
                case 2
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
                otherwise
                    error('Invalid syntax for the "tvumat" command.');
            end
            obj.Data = umat(obj.Data);
            Nt = numel(obj.Time);
            if ~obj.isTimeInvariant
                % XXX-Undocumented: The next code allows Data to be an
                % Ny-by-Nu-by-Nt-by-1 or an Ny-by-Nu-by-1-by-Nt vector
                % where Nt is the dim of Time. It reshapes Data to be Ny-
                % -by-Nu-by-1-by-Nt. This ensures that obj.Data always has
                % at least 4 dims with the last dim corresponding to Time.
                % This simplifies some complexity as Matlab treats an SS
                % array of size Ny-by-Nu-by-Nt-by-1 as having four dims,
                % i.e. it does not drop the trailing "1". As a result,
                % SS arrays have dim of 2 or >=4 but not 3.
                if nmodels(obj.Data)==Nt
                    obj.Data = reshape(obj.Data,[1 Nt]);
                end
            end
            
            % XXX Squeeze extra singleton dimensions?
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
        end
        
        %% ISVALID
        % ISVALID Determine if TVUMAT object is valid.
        function [pflag,errstr] = isvalid(obj)
            errstr = [];
            pflag = 1;
            
            % Check Data Class
            Data = obj.Data; %#ok<*PROP>
            if ~isa(Data,'umat')
                pflag = 0;
                errstr = 'Data must be an umat.';
                return
            end
            
            % Check Time
            Time = obj.Time;
            Nt = numel(Time);
            if ~isa(Time,'double') || Nt<2 || any( diff(Time)< 0 )
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
                error('Discrete-Time TVUMAT are not supported yet.');
            end
            
            % Check Time/Data consistency
            szData = size(obj.Data);
            if obj.isTimeInvariant && ~all(Time==[-inf inf])
                pflag = 0;
                errstr = 'Time vector must be [-inf inf] for "Constant" DataType.';
            end
            if ~obj.isTimeInvariant && ~all(size(Data)==0) && ...
                    ~(numel(szData)>2 && szData(end)==Nt )
                pflag = 0;
                errstr = 'Dimensions of Time and Data are incompatible.';
            end
            if ~isequal(objTs,0) && ~isequal(objTs,-1)
                % Discrete-Time TVUMAT
                diffTime = diff(objTime);
                
                if ~isequal(numel(uniquetol(diffTime,sqrt(eps))),1)
                    pflag = 0;
                    errstr = 'Discrete-Time TVUMAT must be regularly sampled.';
                elseif ~isequal(numel(uniquetol([uniquetol(diffTime,sqrt(eps)),...
                        obj.Ts],sqrt(eps))),1)
                    pflag = 0;
                    errstr = 'SampleTime Property must be consistent with Time.';
                end
            end
            
            % Check IM
            IM = obj.InterpolationMethod;
            if ~isequal(objTs,0)
                % Discrete-Time TVUMAT IM should be empty
                if ~isempty(IM)
                    pflag = 0;
                    errstr = 'InterpolationMethod must be empty for Discrete-Time TVUMAT.';
                end
            elseif ~any(strcmpi(IM,{'Flat'; 'Nearest';'Linear'; 'Spline'}) )
                pflag = 0;
                errstr = ['InterpolationMethod must be "Flat", '...
                    '"Nearest", "Linear", or "Spline" for Continuous-Time TVUMAT'];
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
        
        function out = get.NominalValue(obj)
            out = obj.Data.NominalValue;
        end
        
        function out = get.Uncertainty(obj)
            out = obj.Data.Uncertainty;
        end
        
        function out = get.SamplingGrid(obj)
            out = obj.Data.SamplingGrid;
        end
        
        function out = get.Name(obj)
            out = obj.Data.Name;
        end
        
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
        
        % XXX Implement: numArgumentsFromSubscript
        
        %% SET Methods
        function obj = set.NominalValue(obj,V)
            obj.Data.NominalValue = V;
        end
        
        function obj = set.Uncertainty(obj,V)
            obj.Data.Uncertainty = V;
        end
        
        function obj = set.SamplingGrid(obj,V)
            obj.Data.SamplingGrid = V;
        end
        
        function obj = set.Name(obj,V)
            obj.Data.Name = V;
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
        
        function varargout = evalt(varargin)
            % EVALT Evaluate a TVUMAT(s) at specified time(s)
            %
            % B = EVALT(A,T) evaluates a TVUMAT A on the vector of times T
            % using the interpolation method specified by
            % A.InterpolationMethod.  B is returned as a TVUMAT.
            %
            % [C,D,...] = EVALT(A,B,...,T) evaluates a TVUMAT A and TVUMAT
            % B etc. on the vector of times T. [C,D,...] is returned as
            % TVUMATs.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(true,varargin{:});
        end
        
        function varargout = tvsubs(varargin)
            % TVSUBS Evaluate a TVSUBS(s) at specified time(s)
            %
            % B = TVSUBS(A,T) evaluates a TVUMAT A on the vector of times T
            % using the interpolation method specified by
            % A.InterpolationMethod.  B is returned as a umat array.
            %
            % [C,D,...] = TVSUBS(A,B,...,T) evaluates a TVUMAT A and TVUMAT
            % B etc. on the vector of times T. [C,D,...] is returned as
            % umat arrays.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(false,varargin{:});
        end
        
        function varargout = tvsplit(A,aTf)
            % TVSPLIT splits the TVUMAT on specified time intervals
            %
            % [B,C] = TVSPLIT(A,T1) splits the horizon of the tvumat A and
            % returns B as A in [T0 T1]
            %         C as A in [T1 Tf]
            %
            % [B,C,D,E,...] = TVSPLIT(A,[T1,T2,T3,...]) splits the horizon of the
            % tvumat A and returns B as A in [T0 T1]
            %         C as A in [T1 T2]
            %         D as A in [T2 T3]
            %         E as A in [T3 Tf]
            %
            % NOTE: Here, [T0,Tf] = getHorizon(A) are assumed to be boundry points,
            % input tvumat must be finite horizon and T1, T2, T3 must be sorted such
            % that T0 < T1 < T2 < T3 < ... < Tf.
            %
            % Example:
            % t = 0:0.01:10;
            % A = evalt(tvumat(ureal('A',0)),t);
            % [B,C,D,E] = tvsplit(A,[3 5 8]);
            [varargout{1:max(nargout,1)}] = ltvutil.tvsplit(A,aTf);
        end
        
        function varargout = tvmerge(varargin)
            % TVMERGE merges the TVUMATs on specified horizons
            [varargout{1:max(nargout,1)}] = ltvutil.tvmerge(varargin{:});
        end
        
        function out = length(obj)
            out = size(obj);
            out = max(out);
        end
        
        function out = isempty(obj)
            %szo = size(obj); out = any(~szo);
            out = isempty(obj.Data);
        end
        
        function out = ndims(obj)
            if obj.isTimeInvariant
                out = ndims( obj.Data );
            else
                out = ndims( obj.Data ) - 1;
                szData = size(obj.Data);
                % UMAT Arrays have a min dimension = 4
                if out==3 && szData(end-1)==1
                    % Handles case of Ny-by-Nu-by-1-by-NT
                    out=2;
                elseif out==3
                    % Handles case of Ny-by-Nu-by-Na-by-NT
                    out=4;
                end
            end
        end
        
        function varargout = size(obj,varargin)
            if obj.isTimeInvariant
                [varargout{1:max(nargout,1)}] = size( obj.Data, varargin{:});
            else
                Data = obj.Data; %#ok<*PROPLC>                
                nd = ndims(Data)-3;
                id = repmat({':'},1,nd);
                [varargout{1:max(nargout,1)}] = size( Data(:,:,id{:},1), varargin{:});
            end
        end
        
        function out = isUncertain(obj)
            % Returns true if Uncertain object
            out = isUncertain(obj.Data);
        end
        
        function [M,Del,Blkstruct,NormUnc] = lftdata(obj)
            [M,Del,Blkstruct,NormUnc] = ltvutil.lftdata(obj);
        end
        
        % XXX To implement:
        % iosize, usubs, usample, simplify, connect, wcnorm,
        
        %% Unary Operations: Element-by-Element
        function out = fliplr(obj)
            out = tvunopebe(@flip,2,obj);
        end
        function out = flipud(obj)
            out = tvunopebe(@flip,1,obj);
        end
        function out = uminus(obj)
            out = tvunopebe(@uminus,obj);
        end
        function out = uplus(obj)
            out = tvunopebe(@uplus,obj);
        end
        
        %% Unary Matrix Operations: For-loop over time
        function out = ctranspose(obj)
            out = tvunop(@ctranspose,obj);
        end
        function out = diag(obj,varargin)
            out = tvunop(@diag,obj,varargin{:});
        end
        function out = inv(obj)
            out = tvunop(@inv,obj);
        end
        function out = permute(obj,varargin)
            out = tvunop(@permute,obj,varargin{:});
        end
        function out = repmat(obj,varargin)
            out = tvunop(@repmat,obj,varargin{:});
        end
        function out = reshape(obj,varargin)
            out = tvunop(@reshape,obj,varargin{:});
        end
        function out = transpose(obj)
            out = tvunop(@transpose,obj);
        end
        
        %% Unary Matrix Operations: For-loop + Variable Outputs
        
        
        %% Binary Matrix Operations : Single Output
        function out = feedback(A,B,varargin)
            out = tvbinop(@feedback,A,B,varargin{:});
        end
        function out = lft(A,B,varargin)
            out = tvbinop(@lft,A,B,varargin{:});
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
        function out = plus(A,B)
            out = tvbinop(@plus,A,B);
        end
        
        %% N-ary Matrix Operations
        
        %% Recursive Binary Operations
        function out = append(varargin)
            out = tvrecop(@append,varargin{:});
        end
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
            % Merge with size
            header = ' Uncertain time-varying matrix';
            [Nr,Nc] = size(obj);
            header = [header sprintf(' with %d rows and %d columns.',Nr,Nc)];
            header = [newline header newline newline];
        end
        
        %% GetFooter
        function footer = getFooter(obj)             %#ok<MANU>
            footer = ' Use "get" method to see all the properties.';
            footer = [newline footer newline newline];
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
            % Header
            fprintf('%s',getHeader(obj));
            
            % Print relevant details
            propList = [];
            propList.Data = obj.Data;
            propList.Horizon = [obj.Time(1),obj.Time(end)];
            propList.InterpolationMethod = obj.InterpolationMethod;
            propgrp = matlab.mixin.util.PropertyGroup(propList);
            matlab.mixin.CustomDisplay.displayPropertyGroups(obj,propgrp);
            
            % Footer
            fprintf('%s',getFooter(obj));
        end
    end
end % end of classdef
