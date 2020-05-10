classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvlsimOptions < tvodeOptions
    
    % Options set for TVNORM
    
    properties
        % Specify Type of Simulation
        Type = 'ForwardInTime';
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
        
    end
end
