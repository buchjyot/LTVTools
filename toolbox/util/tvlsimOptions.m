classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvlsimOptions < tvodeOptions
    
    % Options set for TVLSIM
    
    properties
        % Specify Type of Simulation [{'ForwardInTime'} | {'BackwardInTime'}].
        Type = 'ForwardInTime';
        
        % StepSize [{'Default} | {'Auto'} | {0.1}].
        % Default: option uses input time grid as integration time grid
        % Auto   : time grid is determined by ODE solver
        % 0.1    : time grid for integration is T0:0.1:Tf, where T0, Tf is
        % plant horizon
        StepSize = 'Default';
    end
    
    methods
        %% Constructor
        function opt = tvlsimOptions(varargin)
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
        function opt = set.Type(opt,V)
            V = ltipack.matchKey(V,{'ForwardInTime','BackwardInTime'});
            if isempty(V)
                error('The "Type" option must be set to ''ForwardInTime'' or ''BackwardInTime''.')
            end
            opt.Type = V;
        end
        
        %% Specify StepSize
        function opt = set.StepSize(opt,V)
            switch class(V)
                case 'char'
                    V = ltipack.matchKey(V,{'Default','Auto'});
                case 'double'
                    % Do nothing
                otherwise
                    V = [];
            end
            if isempty(V)
                error('The "StepSize" option must be set to ''Default'' or ''Auto'' or double scalar.')
            end
            opt.StepSize = V;
        end
        
    end
end
