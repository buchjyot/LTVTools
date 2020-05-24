classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) slpoweritOptions < tvpoweritOptions
    
    % Options set for SLPOWERIT
    
    properties        
        % Objective
        Objective = 'L2toL2';
        
        % InputL2Norm
        InputL2Norm = 1;
        
        % Linearization Options
        LinOpt = linearizeOptions('BlockReduction','off');
    end
    
    methods
        %% Constructor
        function opt = slpoweritOptions(varargin)
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
       
        %% Specify Objective
        function opt = set.Objective(opt,V)
            V = ltipack.matchKey(V,{'L2toL2','L2toE'});
            if isempty(V)
                error('The "Objective" option must be set to ''L2toL2'' or ''L2toE''.')
            end
            opt.Objective = V;
        end              
        
        %% Specify LinOpt
        function opt = set.LinOpt(opt,V)
            if isa(V,'linearize.LinearizeOptions')
                opt.LinOpt = V;
            else
                error('The "LinOpt" must be set using "linearizeOptions".')
            end
        end
        
        %% Specify InputL2Norm
        function opt = set.InputL2Norm(opt,V)
            if isa(V,'double') && isscalar(V) && V > 0
                opt.InputL2Norm = V;
            else
                error('The "InputL2Norm" must be set to a double scalar.')
            end
        end
        
    end
end