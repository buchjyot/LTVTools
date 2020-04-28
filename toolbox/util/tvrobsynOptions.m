classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvrobsynOptions
    
    properties
        % Display progress
        Display = 'off';
        
        % Tolerances for DK-iteration termination condition
        RelTol = 1e-3;
        AbsTol = 1e-4;
        
        % Flag whether to stop DK-iterations when wcgain within tolerance
        StopWhenWithinTol = false;
        
        % Maximum iterations for DK synthesis
        MaxIter = 10;
        
        % Synthesis and Analysis Options
        SynthesisOptions = tvhinfsynOptions;
        AnalysisOptions = tvwcOptions;
        
        % Initial D-scale
        % Can be a function handle or matrix itself
        InitialDScale = @(Nw,Nv) eye(Nw,Nv);
        
        % Debug Mode
        DebugMode = false;
    end
    
    methods
        %% Constructor
        function opt = tvrobsynOptions(varargin)
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
        
        %% Specify RelTol
        function opt = set.RelTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.RelTol = V;
            else
                error('The "RelTol" option must be a non-negative scalar.')
            end
        end
        
        %% Specify AbsTol
        function opt = set.AbsTol(opt,V)
            if isa(V,'double') && isscalar(V) && V >= 0
                opt.AbsTol = V;
            else
                error('The "AbsTol" option must be a non-negative scalar.')
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
        
        %% Specify SynthesisOptions
        function opt = set.SynthesisOptions(opt,V)
            if isa(V,'tvhinfsynOptions')
                opt.SynthesisOptions = V;
            else
                error('SynthesisOptions must be specified as ''tvhinfsynOptions''.');
            end
        end
        
        %% Specify AnalysisOptions
        function opt = set.AnalysisOptions(opt,V)
            if isa(V,'tvwcOptions')
                opt.AnalysisOptions = V;
            else
                error('AnalysisOptions must be specified through ''tvwcOptions''.');
            end
        end
        
        %% Specify InitialDScale
        function opt = set.InitialDScale(opt,V)
            if isa(V,'function_handle') || isa(V,'double')
                opt.InitialDScale = V;
            else
                error('InitialDScale must be either function_handle or a double.');
            end
        end
        
        %% Set StopDKIterWhenWCGainWithinTol
        function opt = set.StopWhenWithinTol(opt,V)
            opt.StopWhenWithinTol = boolean(V);
        end
        
        %% Set DebugMode
        function opt = set.DebugMode(opt,V)
            opt.DebugMode = boolean(V);
        end
    end
end
