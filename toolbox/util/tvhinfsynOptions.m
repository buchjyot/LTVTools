classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvhinfsynOptions < tvnormOptions
    
    properties
        %% Inherited Properties:
        % Display progress of computation [{'off'} | 'on'].
        % Relative Tolerance for bisection
        % Absolute Tolerance for bisection
        % Initial bounds for bisection
        % OdeSolver Selection
        % OdeOptions Specifications
        
        % Method for synthesis.
        Method = 'MinGamma';
        
        % Multiplicative factor ( >= 1 ) to back off the minimum gamma when
        % Method = 'BackOff'.
        BackOffFactor = 1.2;
    end
    
    methods
        %% Constructor
        function opt = tvhinfsynOptions(varargin)
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
        
        %% Specify Method
        function opt = set.Method(opt,V)
            AllowableVal = {'MinGamma'; 'BackOff';};
            if ischar(V) && any( strcmp(V,AllowableVal) )
                opt.Method = V;
            else
                error('The "Method" must be set to one of the following:\n%s, %s.\n',...
                    AllowableVal{:});
            end
        end
        
        %% Specify BackOffFactor
        function opt = set.BackOffFactor(opt,V)
            if isscalar(V) && isa(V,'double') && V >=1
                opt.BackOffFactor = V;
            else
                error('BackOffFactor must be a scalar >= 1.');
            end
        end
    end
end