classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvwcOptions
    
    % Options set for tvwcgain (which includes LMI step and RDE step)
    
    properties
        % Uncertainty Level
        ULevel = 1;
        
        % RDE and LMI Options
        RDEOptions = tvnormOptions;
        LMIOptions = tvlmiOptions;
        
        % Maximum iterations for analysis
        MaxIter = 10;
        
        % Stopping Tolerance for analysis
        StopTol = 5e-3;
        
        % Display the iteration info
        % NOTE: This is a top level display which will print the iteration
        % info. To print the bisection info refer to RDEOptions display
        Display = 'off';
        
        % Finite Horizon is divided in to Nlmi numbers of linearly spaced
        % time-points for solving LMI on a corse grid
        Nlmi = 20;
        
        % Time-points for Spline basis functions
        Nsp = 10;
    end
    
    methods
        %% Constructor
        function opt = tvwcOptions(varargin)
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
        
        %% Specify ULevel
        function opt = set.ULevel(opt,V)
            % Uncertainty level can be either "double" or "tvmat"
            cond = (V>0); eflag=false;
            switch class(V)
                case 'double'
                    if cond, opt.ULevel = V; else, eflag = true; end
                case 'tvmat'
                    if all(all(cond.Data)), opt.ULevel = V; else, eflag = true; end
            end
            if eflag
                error('The "ULevel" option must be non-negative number or time-varying matrix.');
            end
        end
        
        %% Specify MaxIter
        function opt = set.MaxIter(opt,V)
            if V > 0
                opt.MaxIter = round(V);
            else
                error('The "MaxIter" option must be non-negative integer.');
            end
        end
        
        %% Specify StopTol
        function opt = set.StopTol(opt,V)
            if V > 0
                opt.StopTol = V;
            else
                error('The "StopTol" option must be non-negative value.');
            end
        end
        
        %% Specify Nlmi
        function opt = set.Nlmi(opt,V)
            if isinf(V) || isnan(V)
                error('Number of grid points for LMI must be a natural number.');
            end
            opt.Nlmi = round(V);
        end
        
        %% Specify NSp
        function opt = set.Nsp(opt,V)
            if isinf(V) || isnan(V)
                error('Number of grid points for Spline basis function must be a natural number.');
            end
            opt.Nsp = round(V);
        end
        
        %% Specify RDEOptions
        function opt = set.RDEOptions(opt,V)
            if isa(V,'tvnormOptions')
                opt.RDEOptions = V;
            else
                error('RDEOptions must be specified as ''tvnormOptions''.');
            end
        end
        
        %% Specify LMIOptions
        function opt = set.LMIOptions(opt,V)
            if isa(V,'tvlmiOptions')
                opt.LMIOptions = V;
            else
                error('LMIOptions must be specified as ''tvlmiOptions''.');
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
        
    end
end