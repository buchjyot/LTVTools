classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) poweritOptions < tvodeOptions
    
    % Options set for WCAOPTIONS
    
    properties
        % Display progress of computation [{'on'} | {'plot'} | {'off'}].
        Display = 'off';
        
        % Stopping Tolerance
        StopTol = 1e-2;
        
        % Maximum Iterations
        MaxIter = 500;                
        
        % PartialsFunc
        PartialsFunc = [];
        
        % StoreAllIterations
        % By default this is false because, sometimes this iteration info
        % can require a lot of memory. We do not store all the information
        % unless user wants it.
        StoreAllIter = false;
        
        % Simulation Options
        SimOpt = simset('RelTol',1e-3,'AbsTol',1e-6);
        
        % Linearization Options
        LinOpt = linearizeOptions('BlockReduction','off');
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
            V = ltipack.matchKey(V,{'off','plot','on'});
            if isempty(V)
                error('The "Display" option must be set to ''on'', ''off'' or ''plot''.')
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
        
        %% Specify LinOpt
        function opt = set.LinOpt(opt,V)
            if isa(V,'linearize.LinearizeOptions')
                opt.LinOpt = V;
            else
                error('The "LinOpt" option must be set using ''linearizeOptions''.')
            end
        end
        
        %% Specify SimOpt
        function opt = set.SimOpt(opt,V)
            if isa(V,'struct')
                opt.SimOpt = V;
            else
                error('The "SimOpt" option must be a structure specified using ''simset''.')
            end
        end
        
        %% Specify PartialsFunc
        % PartialsFunc  Function that returns the partial derivatives of
        % the system function f(x,u,t) and the output function g(x,u,t).
        % The calling syntax is [dfdx,dfdu,dgdx,dgdu] =
        % function_name(x,u,t).
        function opt = set.PartialsFunc(opt,V)
            if isa(V,'function_handle')
                opt.PartialsFunc = V;
            else
                error('The "PartialsFunc" option must be set to a function_handle');
            end
        end
        
        %% ReturnAllIter
        % Specify true if you want to return all the iterations info in
        % the third output argument of powerit
        function opt = set.StoreAllIter(opt,V)
            opt.StoreAllIter = boolean(V);
        end
    end
end