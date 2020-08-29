%% TVSS create a time-varying state-space system
%
%   G = tvss(Data,Time) creates a time-varying system. Time is an Nt-by-1
%   vector of increasing time values. If Data is a 1-by-Nt SS array then G
%   is a single-input, single-output time-varying system with value Data(i)
%   at Time(i). If Data is N1-by-N2-by-Nt then G is an N1-by-N2 time-
%   varying system with value Data(:,:,i) at Time(i). In general, if Data
%   is an N1-by-N2-by-...-by-Nad-by-Nt array then G is an N1-by-N2-by-...
%   -by-Nad time-varying SS array with value Data(:,...,:,i) at Time(i).
%
%   G = tvss(Data,Time,InterpolationMethod) specifies the method to
%   interpolate the state-space matrices between the specified time grid
%   points.  The method can be 'Linear', 'Flat', 'Nearest', or 'Spline'.
%   Default is InterpolationMethod = 'Linear'.
%
%   G = tvss(Data) generates a constant system.
%
%   G = tvss(A,B,C,D) creates a continuous-time, time-varying system with
%   state matrices defined by TVMATs (A,B,C,D). Additional model properties
%   can be set using name/value pairs.  See help for SS command.
%
%   G = tvss(A,B,C,D,Ts) creates a discrete-time, time-varying system with
%   sample time Ts.
%
%   % EXAMPLE: (CUT/PASTE)
%   % Create 1-state system with time-varying "A" matrix.
%   Time = linspace(0,5,10)';
%   AData = -5+0.1*Time.^2;
%   A = tvmat(AData,Time);
%   B = 1; C = 1; D=0;
%   G = tvss(A,B,C,D);
%
% See also: tvmat

classdef (InferiorClasses={?frd, ?ss,?tf,?zpk,?ureal,?ucomplex,?ucomplexm,...
        ?ultidyn,?udyn,?umat,?uss,?ufrd,?tvmat}) tvss < matlab.mixin.CustomDisplay
    % Class definition for time-varying dynamical system
    
    properties (SetAccess = public)
        Data = ss;
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
        % Following properties are dependent on underlying ssarray object
        % High-Level Summary of Low-Level SSARRAY Data
        Ts
        TimeUnit
        Name
        StateName
        StateUnit
        InputName
        InputUnit
        InputGroup
        OutputName
        OutputUnit
        OutputGroup
        Notes
        UserData
    end
    
    methods
        %% Constructor
        % XXX Add undocumented flag to skip error checking?
        % XXX Currently disallows descriptor data.
        function obj = tvss(Data,Time,InterpolationMethod,varargin)
            nin = nargin;
            switch nin
                case 0
                    % Return default object
                    return
                case 1
                    if isa(Data,'tvss')
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
                otherwise
                    % tvss(A,B,C,D,Ts,...) where (A,B,C,D) are TVMATs or
                    % TVUMATs
                    [isTimeInvariant,A,B,C,D] = tvswitchyard(Data,Time,InterpolationMethod,varargin{1});
                    varargin = varargin(2:end);
                    
                    if isTimeInvariant
                        obj.Data = ss(A.Data,B.Data,C.Data,D.Data,varargin{:});
                        obj.isTimeInvariant = true;
                    else
                        obj.Data = ss(A.Data,B.Data,C.Data,D.Data,varargin{:});
                        obj.Time = A.Time;
                        obj.InterpolationMethod = A.InterpolationMethod;
                    end
            end
            
            % Convert Data in to State-Space
            obj.Data = ss(obj.Data);
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
                elseif nmodels(obj.Data)==1
                    % Expand across the time dimention
                    obj.Data = repsys(obj.Data,[1 1 1 Nt]);
                end
            end
            
            % Use isvalid to perform all error checking
            [pflag,errstr] = isvalid(obj);
            if pflag==0
                error(errstr);
            end
            % XXX Squeeze extra singleton dimensions?
        end
        
        %% ISVALID
        function [pflag,errstr] = isvalid(obj)
            % ISVALID Determine if TVSS object is valid.
            errstr = [];
            pflag = 1;
            
            % Check Data Class
            objData = obj.Data;
            if ~isa(objData,'ss')
                pflag = 0;
                errstr = 'Data must be an ss.';
                return
            end
            
            % Check Consistent State Order
            if ~all(order(objData)==order(obj))
                pflag = 0;
                errstr = 'State order must be constant across the time grid.';
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
            
            % Check Time/Data consistency
            szData = size(objData);
            if any(isnan(objTime))
                pflag = 0;
                errstr = 'Time vector must not contain any NaNs.';
            end
            if obj.isTimeInvariant && ~all(objTime==[-inf inf])
                pflag = 0;
                errstr = 'Time vector must be [-inf inf] for "Constant" DataType.';
            end
            if ~obj.isTimeInvariant && ~all(size(objData)==0) && ...
                    ~(numel(szData)>2 && szData(end)==Nt )
                pflag = 0;
                errstr = 'Dimensions of Time and Data are incompatible.';
            end
            
            % Check IM
            IM = obj.InterpolationMethod;
            objTs = obj.Ts;
            if ~isequal(objTs,0)
                error('Discrete-Time TVSS are not supported yet.');
            end
            if ~isequal(objTs,0)
                % Discrete-Time TVSS IM should be empty
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
        
        %% GET Methods
        function V = get.Ts(obj)
            V = obj.Data.Ts;
        end
        
        function V = get.TimeUnit(obj)
            V = obj.Data.TimeUnit;
        end
        
        function V = get.StateName(obj)
            V = unique(obj.Data.StateName);
        end
        
        function V = get.StateUnit(obj)
            V = unique(obj.Data.StateUnit);
        end
        
        function V = get.InputName(obj)
            V = unique(obj.Data.InputName);
        end
        
        function V = get.InputUnit(obj)
            V = unique(obj.Data.InputUnit);
        end
        
        function V = get.InputGroup(obj)
            V = obj.Data.InputGroup;
        end
        
        function V = get.OutputName(obj)
            V = unique(obj.Data.OutputName);
        end
        
        function V = get.OutputUnit(obj)
            V = unique(obj.Data.OutputUnit);
        end
        
        function V = get.OutputGroup(obj)
            V = obj.Data.OutputGroup;
        end
        
        function V = get.Notes(obj)
            V = obj.Data.Notes;
        end
        
        function V = get.UserData(obj)
            V = obj.Data.UserData;
        end
        
        function V = get.Name(obj)
            V = obj.Data.Name;
        end
        
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
            obj.Data.Ts = V;
        end
        
        function obj = set.TimeUnit(obj,V)
            obj.Data.TimeUnit = V;
        end
        
        function obj = set.StateName(obj,V)
            obj.Data.StateName = V;
        end
        
        function obj = set.StateUnit(obj,V)
            obj.Data.StateUnit = V;
        end
        
        function obj = set.InputName(obj,V)
            obj.Data.InputName = V;
        end
        
        function obj = set.InputUnit(obj,V)
            obj.Data.InputUnit = V;
        end
        
        function obj = set.InputGroup(obj,V)
            obj.Data.InputGroup = V;
        end
        
        function obj = set.OutputName(obj,V)
            obj.Data.OutputName = V;
        end
        
        function obj = set.OutputUnit(obj,V)
            obj.Data.OutputUnit = V;
        end
        
        function obj = set.OutputGroup(obj,V)
            obj.Data.OutputGroup = V;
        end
        
        function obj = set.Notes(obj,V)
            obj.Data.Notes = V;
        end
        
        function obj = set.UserData(obj,V)
            obj.Data.UserData = V;
        end
        
        function obj = set.Name(obj,V)
            obj.Data.Name = V;
        end
        
        function out = set(obj,varargin)
            % SET  Modifies values of dynamic system properties.
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
        % connect, getLoopTransfer, getIOTransfer, numArgumentsFromSubscript
        
        function out = ss(obj)
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
        
        function out = isUncertain(obj)
            % Returns true if Uncertain object
            out = isUncertain(obj.Data);
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
        
        function out = isempty(obj)
            out = any(size(obj)==0);
        end
        
        function out = ndims(obj)
            if obj.isTimeInvariant
                out = ndims( obj.Data );
            else
                out = ndims( obj.Data ) - 1;
                szData = size(obj.Data);
                % SS Arrays have a min dimension = 4
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
            nout = nargout;
            if obj.isTimeInvariant
                [varargout{1:max(nout,1)}] = size( obj.Data, varargin{:});
            else
                Data = obj.Data; %#ok<*PROPLC>
                
                nd = ndims(Data)-3;
                id = repmat({':'},1,nd);
                [varargout{1:max(nout,1)}] = size( Data(:,:,id{:},1), varargin{:});
            end
        end
        
        function [GT,Ti,Tdot] = ss2ss(G,T)
            % SS2SS  Change of state coordinates for TVSS.
            %
            %   SYS = SS2SS(SYS,T) performs the similarity transformation
            %   z = Tx on the state vector x of the time-varying state-space
            %   model SYS. The resulting state-space model is described by:
            %      .               -1
            %      z = [Tdot+T*A]*T   z + TB u
            %                      -1
            %               y = CT   z + D u
            
            % Differentiate TVMAT (Force to Spline InterpolationMethod)
            Td = T;
            Td.InterpolationMethod = 'Spline';
            Tdot = tvdiff(Td);
            Tdot.InterpolationMethod = G.InterpolationMethod;
            
            % Similarity Transformation
            Ti = 1/T;
            [A,B,C,D] = ssdata(G);
            AT = (Tdot + T*A)*Ti;
            BT = T*B;
            CT = C*Ti;
            DT = D;
            GT = tvss(AT,BT,CT,DT);
        end
        
        function out = order(obj)
            out = order(obj.Data(:,:,1));
        end
        
        function out = eig(obj)
            A = ssdata(obj);
            out = eig(A);
        end
        
        function varargout = evalt(varargin)
            % EVALT Evaluate a TVSS on a specified time grid
            %
            % H = EVALT(G,T) evaluates a TVSS G on the vector of times T
            % using the interpolation method specified by
            % G.InterpolationMethod.  H is returned as a TVSS object.
            %
            % [GT,KT,...] = EVALT(G,K,...,T) evaluates a TVSS G and K etc.
            % on the vector of times T. GT,KT,... are returned as TVSS
            % objects.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(true,varargin{:});
        end
        
        function varargout = tvsubs(varargin)
            % TVSUBS Evaluate a TVSS on a specified time grid
            %
            % H = TVSUBS(G,T) evaluates a TVSS G on the vector of times T
            % using the interpolation method specified by
            % G.InterpolationMethod.  H is returned as ss array.
            %
            % [GT,KT,...] = TVSUBS(G,K,...,T) evaluates a TVSS G and K etc.
            % on the vector of times T. GT,KT,... are returned as ss array.
            [varargout{1:max(nargout,1)}] = ltvutil.evalt(false,varargin{:});
        end
        
        function varargout = tvsplit(A,aTf)
            % TVSPLIT splits the TVSS on specified time intervals
            %
            % [B,C] = TVSPLIT(A,T1) splits the horizon of the tvss A and
            % returns B as A in [T0 T1]
            %         C as A in [T1 Tf]
            %
            % [B,C,D,E,...] = TVSPLIT(A,[T1,T2,T3,...]) splits the horizon
            % of the tvss A and returns B as A in [T0 T1]
            %         C as A in [T1 T2]
            %         D as A in [T2 T3]
            %         E as A in [T3 Tf]
            %
            % NOTE: Here, [T0,Tf] = getHorizon(A) are assumed to be boundry
            % points, input tvss must be finite horizon and T1, T2, T3 must
            % be sorted such that T0 < T1 < T2 < T3 < ... < Tf.
            %
            % Example:
            % t = 0:0.01:10;
            % A = evalt(tvss(rss(4,3,2)),t);
            % [B,C,D,E] = tvsplit(A,[3 5 8]);
            [varargout{1:max(nargout,1)}] = ltvutil.tvsplit(A,aTf);
        end
        
        function varargout = tvmerge(varargin)
            % TVMERGE merges the TVSS on specified horizons
            [varargout{1:max(nargout,1)}] = ltvutil.tvmerge(varargin{:});
        end
        
        %% Synthesis functions
        function [K,CL,GAM,INFO] = tvsyn(G,Ny,Nu,varargin)
            % Time-Varying Hinfinity Synthesis
            % If Ny == 0 then performs state-feedback synthesis
            % If Ny ~= 0 then performs measurement-feedback synthesis
            
            % Minimum 3 and maximum 6 input arguments are allowed
            narginchk(3,6);
            
            switch Ny
                case 0
                    % NOTE: this assumes that user provided plant G suitable for full
                    % state feedback synthesis, try >> help tvhinfsfb for details
                    [K,CL,GAM,INFO] = tvhinfsfb(G,Nu,varargin{:});
                otherwise
                    [K,CL,GAM,INFO] = tvhinfsyn(G,Ny,Nu,varargin{:});
            end
        end
        
        function [Nx,Ne,Nd,NL2,NY,NU] = syniodim(G,Ny,Nu,NE)
            % Returns IO and state dimentions for synthesis in the
            % specific order
            if nargin < 4, NE = 0; end
            [NY,NU] = size(G);
            Nx = order(G);
            Ne = NY - Ny;
            NL2 = Ne - NE;
            Nd = NU - Nu;
        end
        
        %% Unary Operations: Element-by-Element
        function out = conj(obj)
            out = tvunopebe(@conj,obj);
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
        function out = inv(obj)
            out = tvunop(@inv,obj);
        end
        function out = norm(obj,varargin)
            out = tvunop(@norm,obj,varargin{:});
        end
        function out = permute(obj,varargin)
            out = tvunop(@permute,obj,varargin{:});
        end
        function out = repmat(obj,varargin)
            out = tvunop(@repmat,obj,varargin{:});
        end
        function out = repsys(obj,varargin)
            out = tvunop(@repsys,obj,varargin{:});
        end
        function out = reshape(obj,varargin)
            out = tvunop(@reshape,obj,varargin{:});
        end
        function out = transpose(obj)
            out = tvunop(@transpose,obj);
        end
        
        %% Unary Matrix Operations: For-loop + Variable Outputs
        
        %% Binary Matrix Operations: Single Output
        function out = feedback(A,B,varargin)
            out = tvbinop(@feedback,A,B,varargin{:});
        end
        function out = ldivide(A,B)
            out = tvbinop(@ldivide,A,B);
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
        function out = rdivide(A,B)
            out = tvbinop(@rdivide,A,B);
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
            % Identify Type of TVMAT
            if isequal(obj.Ts,0)
                header = ' Continuous-time time-varying state-space model';
            elseif ~isequal(obj.Ts,0) && ~isequal(obj.Ts,-1) && ~isempty(obj)
                header = ' Discrete-time time-varying state-space model';
            end
            
            % Merge with size
            [Nr,Nc] = size(obj);
            Nx = order(obj);
            header = [header sprintf(' with %d outputs, %d inputs and %d states.',Nr,Nc,Nx)];
            
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
            [sz1,sz2] = size(obj);
            if isequal(sz1,0) && isequal(sz2,0)
                fprintf(newline);
                fprintf(' %s','Empty time-varying state-space model.');
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