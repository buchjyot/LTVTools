classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvodeOptions
    
    % ODE Options set for time varying ODEs
    
    % Example:
    % >> cdreOpt = tvodeOptions;
    % >> cdreOpt.OdeOptions = odeset('RelTol',1e-5,'AbsTol',1e-8);
    
    properties
        % NOTE: Following are also the default for Simulink
        % Specify default solver
        OdeSolver = 'ode45';
        % Options for ODE solver
        OdeOptions = odeset('RelTol',1e-3,'AbsTol',1e-6);
    end
    
    methods
        
        %% Constructor
        function opt = tvodeOptions(varargin)
            % Constructor
            narginchk(0,4);
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
        
        %% Specify ODE Solver
        function opt = set.OdeSolver(opt,V)
            % Update this list if new solve becomes available
            availableODESolvers = ...
                {'ode45','ode23','ode113','ode23t','ode23s',...
                'ode23tb','ode15s','ode15i'};
            if any(strcmpi(availableODESolvers, V))
                opt.OdeSolver = lower(V);
            else
                error('ltvtools:tvodeOptions:OdeSolver',...
                    'The "OdeSolver" must be set to one of the following:\n%s, %s, %s, %s, %s, %s, %s, %s.\n',...
                    availableODESolvers{:});
            end
        end
        
        %% Specify ODE Options
        function opt = set.OdeOptions(opt,V)
            if isstruct(V) && any(isfield(V,fieldnames(odeset)))
                opt.OdeOptions = V;
            else
                error('The "OdeOptions" must be a structure with fields specified by odeset.');
            end
        end
        
    end
end