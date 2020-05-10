classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) poweritOptions < tvodeOptions
    
    % Options set for POWERIT
    
    properties
        % Display progress of computation [{'off'} | {'on'}].
        Display = 'off';
        % Stopping Tolerance
        StopTol = 1e-3;
        % Maximum Iterations
        MaxIter = 100;
        % Grid Points
        Ngrid = 500;
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
                opt.RelTol = V;
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
        
        %% Specify Bounds
        function opt = set.Ngrid(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.Ngrid = round(V);
            else
                error('The "Ngrid" option must be a non-negative scalar.')
            end
        end
    end
end