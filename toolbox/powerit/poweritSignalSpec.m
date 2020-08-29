classdef (CaseInsensitiveProperties = true,TruncatedProperties = true) poweritSignalSpec
    
    % This option specify signal specifications for power iterations
    
    properties
        % Number of Euclidean outputs
        NE = 0;
        
        % Number of uncertain output for dynamic norm-bounded uncertainty
        Nv = 0;
        
        % Number of uncertain plant inputs for dynamic norm-bounded uncertainty
        Nw = 0;
        
        % Number of uncertain parameters
        Np = 0;
        
        % Initial input
        InitialInput = 'randn';
        
        % Input signal norm bound
        InputL2Norm = 1;
        
        % Initial conditions
        InitialConditions = 'fixed';
        
        % Initial condition cost matrix
        % Small value of this means, larger uncertainty
        % Bigger value means, we are certain and InitialConditions are
        % known completly. i.e. InitialConditions = 'fixed';
        InitialCondCostMat = @(Nx) diag(inf(Nx,1));
        
        % Dynamic uncertainty norm bound
        UncNormBnd = 1;
        
        % Range of uncertain parameter (Np by 2) vector
        UncParamRange = [0 0];
    end
    
    methods
        %% Constructor
        function opt = poweritSignalSpec(varargin)
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
        
        %% Specify NE
        function opt = set.NE(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.NE = round(V);
            else
                error('The "NE" must be set to a double scalar.')
            end
        end
        
        %% Specify Nv
        function opt = set.Nv(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.Nv = round(V);
            else
                error('The "Nv" must be set to a double scalar.')
            end
        end
        
        %% Specify Nw
        function opt = set.Nw(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.Nw = round(V);
            else
                error('The "Nw" must be set to a double scalar.')
            end
        end
        
        %% Specify Np
        function opt = set.Np(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.Np = round(V);
            else
                error('The "NE" must be set to a double scalar.')
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
        
        %% Specify UncParamRange
        function opt = set.UncParamRange(opt,V)
            if isa(V,'double')
                opt.UncParamRange = V;
            else
                error('The "UncParamRange" must be set to a double (Np by 2) vector.')
            end
        end
        
        %% Specify UncNormBnd
        function opt = set.UncNormBnd(opt,V)
            if isa(V,'double') && isscalar(V) && V > 0
                opt.UncNormBnd = V;
            else
                error('The "UncNormBnd" must be set to a double scalar.')
            end
        end
        
        %% Specify InitialInput
        function opt = set.InitialInput(opt,V)
            V = ltipack.matchKey(V,{'randn','ones','rand'});
            if isempty(V)
                error('The "InitialInput" option must be set to ''randn'', ''rand'' or ''ones''.')
            end
            opt.InitialInput = V;
        end
        
        %% Specify InitialConditions
        function opt = set.InitialConditions(opt,V)
            V = ltipack.matchKey(V,{'fixed','free'});
            if isempty(V)
                error('The "InitialConditions" option must be set to ''fixed'' or ''free''.')
            end
            opt.InitialConditions = V;
        end
        
        %% Specify InitialCondCostMat
        function opt = set.InitialCondCostMat(opt,V)
            if isa(V,'functional_handle')
                opt.InitialCondCostMat = V;
            elseif isa(V,'double')
                [nr,nc] = size(V);
                if isequal(nr,nc)
                    opt.InitialCondCostMat = V;
                else
                    error('The "InitialCondCostMat" option must be set to double square matrix or function_handle.');
                end
            else
                error('The "InitialCondCostMat" option must be set to double square matrix or function_handle.');
            end
        end
    end
end