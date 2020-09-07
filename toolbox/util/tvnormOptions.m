classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvnormOptions < tvodeOptions
    
    % Options set for TVNORM
    
    properties
        % Display progress of computation [{'off'} | {'on'}].
        Display = 'off';
        % Relative Tolerance
        RelTol = 1e-2;
        % Absolute Tolerance
        AbsTol = 1e-4;
        % Initial bounds for bisection
        Bounds = [0, inf];
    end
    
    methods
        %% Constructor
        function opt = tvnormOptions(varargin)
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
        
        %% Specify RelTol
        function opt = set.RelTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.RelTol = V;
            else
                error('The "RelTol" option must be a non-negative scalar.')
            end
        end
        
        %% Specify AbsTol
        function opt = set.AbsTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.AbsTol = V;
            else
                error('The "AbsTol" option must be a non-negative scalar.')
            end
        end
        
        %% Specify Bounds
        function opt = set.Bounds(opt,V)
            if isa(V,'double') && numel(V)==2 && V(2) >= V(1) && ...
                    V(1)>=0 && isfinite(V(1))
                opt.Bounds = V;
            else
                error(['The "Bounds" option must be 1-by-2 vector' ...
                    ' with Bounds(2)>=Bounds(1)>=0'])
            end
        end
    end
end
