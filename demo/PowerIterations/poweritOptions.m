classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) poweritOptions < tvodeOptions
    
    % Options set for POWERIT
    
    properties
        % Display progress of computation [{'off'} | {'on'}].
        Display = 'off';
        
        % Stopping Tolerance
        StopTol = 1e-2;
        
        % Maximum Iterations
        MaxIter = 500;
        
        % StepSize
        StepSize = 'Auto';
        
        % Initial Input
        InitialInput = 'rand';
        
        % Objective
        Objective = 'L2toL2';
    end
    
    methods
        %% Constructor
        function opt = poweritOptions(varargin)
            % Constructor
            narginchk(0,numel(properties(opt))*2);
            nin = nargin;
            
            if isequal(nin,0)
                return
            elseif isequal(mod(nin,2),0)
                % Set Name value pairs
                for i = 1:2:nin
                    opt.(varargin{i}) = varargin{i+1};
                end
            else
                error('Input must be in the form of name-value pairs.');
            end
        end
        
        %% Specify Display
        function opt = set.Display(opt,V)
            V = ltipack.matchKey(V,{'off','on'});
            if isempty(V)
                error('The "Display" option must be set to ''on'' or ''off''.')
            end
            opt.Display = V;
        end
        
        %% Specify StopTol
        function opt = set.StopTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.StopTol = V;
            else
                error('The "StopTol" option must be a non-negative scalar.')
            end
        end
        
        %% Specify MaxIter
        function opt = set.MaxIter(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.MaxIter = round(V);
            else
                error('The "MaxIter" option must be a non-negative scalar.')
            end
        end
        
        %% Specify StepSize
        function opt = set.StepSize(opt,V)
            switch class(V)
                case 'char'
                    V = ltipack.matchKey(V,{'Auto'});
                case 'double'
                    % Do nothing
                otherwise
                    V = [];
            end
            if isempty(V)
                error('The "StepSize" option must be set to ''Auto'' or double scalar.')
            end
            opt.StepSize = V;
        end
        
        %% Specify Initial Input
        function opt = set.InitialInput(opt,V)
            switch class(V)
                case 'char'
                    V = ltipack.matchKey(V,{'ones','rand','randn'});
                case 'tvmat'
                    % Do nothing
                otherwise
                    V = [];
            end
            if isempty(V)
                error('The "InitialInput" option must be set to ''ones'',''rand'',''randn'' or a time-varying matrix specified using tvmat.')
            end
            opt.InitialInput = V;
        end
        
        %% Specify Objective
        function opt = set.Objective(opt,V)
            V = ltipack.matchKey(V,{'L2toL2','L2toE'});
            if isempty(V)
                error('The "Objective" option must be set to ''L2toL2'' or ''L2toE''.')
            end
            opt.Objective = V;
        end
        
        %% Validate Initial Input
        % This function is used to construct or validate the initial input required
        % for the power iteration algorithm.
        function out = validateInitialInput(Opt,Tgrid,Nu)
            I = Opt.InitialInput;
            switch class(I)
                case 'tvmat'
                    % If Input is already TVMAT then verify horizon and input dimention
                    [T0,Tf] = getHorizon(I);
                    [Nr,Nc] = size(I);
                    if ~isequal(Nr,Nu) || ~isequal(Nc,1)
                        error('Initial input dimentions must agree with input dimentions of the model.');
                    end
                    if ~isequal(T0,Tgrid(1)) || ~isequal(Tf,Tgrid(end))
                        error('Initial input must be defined on the same horizon as plant.')
                    end
                    out = I/tvnorm(I);
                    
                case 'char'
                    % Construct input that is either "rand", "randn" or "ones"
                    fh  = str2func(I);
                    Nt  = length(Tgrid);
                    U = tvmat(fh(Nu,1,Nt),Tgrid);
                    out = U/tvnorm(U);
            end
        end
        
    end
end