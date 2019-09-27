classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) tvlmiOptions
    
    % LMI Options set for LMIs
    
    % Example:
    % >> lmiOpt = tvlmiOptions;
    % >> lmiOpt.LmiSolver = 'lmilab';
    
    properties
        % Specify default solver
        LmiSolver = 'lmilab';
        
        % Option for LMI solver LmiOptions gives access to certain control
        % parameters of the optimization code. In mincx, this is a
        % five-entry vector organized as follows:
        %
        % options(1) sets the desired relative accuracy on the optimal
        % value lopt (default = 10–2).
        %
        % options(2) sets the maximum number of iterations allowed to be
        % performed by the optimization procedure (100 by default).
        %
        % options(3) sets the feasibility radius. Its purpose and usage are
        % as for feasp.
        %
        % options(4) helps speed up termination. If set to an integer value
        % J > 0, the code terminates when the objective cTx has not
        % decreased by more than the desired relative accuracy during the
        % last J iterations.
        %
        % options(5) = 1 turns off the trace of execution of the
        % optimization procedure. Resetting options(5) to zero (default
        % value) turns it back on.
        %
        % Setting option(i) to zero is equivalent to setting the
        % corresponding control parameter to its default value. See feasp
        % for more detail.
        
        % Specify default options
        % opts = zeros(5,1);
        % opts(1) = 1e-4;% Relative Accuracy
        % opts(2) = 200; % Max # of iters
        % opts(5) = 1;   % Toggle display
        LmiSolverOptions = [0 200 0 0 1]';
        
        % Initial decision variables for LMI solver
        LmiSolverInit = [];
    end
    
    methods
        
        %% Constructor
        function opt = tvlmiOptions(varargin)
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
        
        %% Specify LMI Solver
        function opt = set.LmiSolver(opt,V)
            % Update this list if new solve becomes available
            AllowableVal = {'lmilab'};
            % XXX: Possible Future Support:
            % 'sedumi','sdpam'; 'csdp'; 'dsdp';'sdpt3';'sdplr';
            if any(strcmpi(AllowableVal, V))
                opt.LmiSolver = lower(V);
            else
                % error('The "LmiSolver" must be set to one of the following:\n%s.\n',...
                % AllowableVal{:});
                error('Only ''lmilab'' is currently supported as "LmiSolver".');
            end
        end
        
        %% Specify Solver Specific LMI Options
        function opt = set.LmiSolverOptions(opt,V)
            opt.LmiSolverOptions = V;
        end
        
        %% Specify LMI Solver Initial Decision Variables
        function opt = set.LmiSolverInit(opt,V)
            opt.LmiSolverInit = V(:);
        end
    end
end