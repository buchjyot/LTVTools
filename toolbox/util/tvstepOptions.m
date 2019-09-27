classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvstepOptions < tvodeOptions & ltioptions.step
    
    properties
        StepTime = 0;
    end
    
    methods
        %% Constructor
        function opt = tvstepOptions(varargin)
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
        
        %% Specify StepTime
        function opt = set.StepTime(opt,V)
            if isa(V,'double') && isscalar(V)
                opt.StepTime = V;
            else
                error('The "StepTime" property must be scalar double.');
            end
        end
    end
end