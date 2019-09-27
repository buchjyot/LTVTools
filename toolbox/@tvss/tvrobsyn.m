function [K,CL,g,info] = tvrobsyn(Gunc,Ny,Nu,varargin)
%% TVROBSYN  Robust controller design on finite horizon using DK-iterations
%
%   [K,CL,GAM,DKINFO] = TVROBSYN(P,NMEAS,NCONT,NEUCLIDEAN,OPT) synthesizes a
%   controller K for the open-loop plant model P via the D-K or D-G-K
%   iteration approach to LTV-synthesis.  P is an uncertain time-varying
%   state-space model (see TVUSS) and it is assumed that the last NMEAS
%   outputs and NCONT inputs of P are the measurement and control channels,
%   respectively. NEUCLIDEAN is the outputs to be penalilzed in a Euclidean
%   Sense.  TVROBSYN returns the controller K, the closed-loop model CL, and
%   the robust closed-loop performance gain GAM These variables are related
%   by:
%
%      CL = lft(P,K);
%      g  = tvwcgain(CL,NE,OPT);

%% Initialize & Input Processing
NE = 0;
narginchk(3,5);
nin = nargin;
% nout = nargout;
switch nin
    case 4
        NE = varargin{1};
        Opts = tvrobsynOptions;
    case 5
        % Default Order
        NE = varargin{1};
        Opts = varargin{2};
end

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

% XXX Check that Euclidean penalty is valid

%% TVROBSYN Engine
% RDE based Synthesis and RDE/DLMI based Analysis
nout = nargout;

% Set Special Flags
DEBUG = true;
AbortFlag = false;

% Read Options
NIterMax = Opts.MaxIter;
DispFlag = isequal(Opts.Display,'on');
tvhopt = Opts.SynthesisOptions;
tvwcopt = Opts.AnalysisOptions;
if DEBUG
    tvnopt = tvnormOptions('Display',tvhopt.Display,'RelTol',tvhopt.RelTol,...
        'AbsTol',tvhopt.AbsTol,'Bounds',tvhopt.Bounds,'OdeSolver',tvhopt.OdeSolver,...
        'OdeOptions',tvhopt.OdeOptions); %#ok<*UNRCH>
end
TermDKIterFlag = Opts.StopWhenWithinTol;
Count = 0;
if TermDKIterFlag
    RelTol = Opts.RelTol;
    AbsTol = Opts.AbsTol;
end
IQCParam = Gunc.UserData;

% Create an extended plant
switch Ny
    case 0
        % Generate Extended System based off Gunc
        [Nx,~,Nd,~] = syniodim(Gunc,0,Nu,NE);
        [Au,Bu,Cu,Du] = ssdata(Gunc);
        Cext = [Cu;eye(Nx)];
        Dext = [Du;zeros(Nx,Nd+Nu)];
        Gext = tvss(Au,Bu,Cext,Dext);
    otherwise
        Gext = Gunc;
end

%% Memory Allocation
Ks       = cell(NIterMax,1);
CLs      = cell(NIterMax,1);
CLoop    = cell(NIterMax,1);
gs       = zeros(NIterMax,1);
wcgain   = zeros(NIterMax,1);
d        = cell(NIterMax,1);
Gscl     = cell(NIterMax,1);
skipinfo = (nout~=4);
if ~skipinfo
    wcinfo  = cell(NIterMax,1);
    infoK   = cell(NIterMax,1);
    info    = cell(NIterMax,1);
end
if DEBUG
    tvnCLs  = cell(NIterMax,1);
end

%% Initialization
% Uncertain Dimentions (Currently only SISO Delta supported)
% XXX This will be part of tvuss later
Nw = 1;
Nv = 1;

% User can specify either function handle or double value for this
Dinit = Opts.InitialDScale;
if isa(Dinit,'function_handle')
    Dinit = Dinit(Nw,Nv);
end

% Gscl stands for scalled open-loop plant
[NY,NU] = size(Gunc); d{1} = Dinit;
Gscl{1} = blkdiag(Dinit,eye(NY-Nv))*Gunc*blkdiag(1/Dinit,eye(NU-Nw));

% Start timing
t0 = tic;

try
    
    %% DK-Iterations
    for i = 1:NIterMax
        
        % Start Display
        if DispFlag
            fprintf(newline);
            fprintf(' ===> DK Iteration #%d:\n',i);
            fprintf(' ### Synthesis Step: \n');
        end
        
        % K-Step
        [Ks{i},CLs{i},gs(i),infoK{i}] = tvsyn(Gscl{i},Ny,Nu,NE,[],tvhopt);
        if DispFlag, fprintf(' Synthesis Gain = %.4f \n',gs(i)); end
        
        % For Debugging K-step verification
        if DEBUG
            gtvn = tvnorm(CLs{i},NE,tvnopt);
            fprintf(' Synthesis TVN = %.4f\n',gtvn(2));
        end
        
        % Plot Controller Gains
        % if DEBUG && isequal(Ny,0)
        % figure(1);clf;tvplot(-Ks{i},'b','LineWidth',3);drawnow;
        % end
        
        % Formulate Closed Loop on Original "Gext" System
        CLoop{i} = lft(evalt(Gext,Ks{i}.Time),Ks{i});
        CLoop{i}.UserData = IQCParam;
        
        % D-Step
        if DispFlag
            fprintf(' ### IQC Analysis Step: \n');
        end
        [wcgain(i),wcinfo{i}] = tvwcgain(CLoop{i},NE,tvwcopt);
        if DispFlag
            fprintf(' IQC Analysis Gain: = %.4f \n', wcgain(i));
        end
        
        % D*D is equal to Psiv'*X11*Psiv
        [~,idx] = min(cellfun(@(in) in.RDE.Upper.Gain,wcinfo{i}));
        X11 = wcinfo{i}{idx}.LMI.X11;
        d{i+1} = sqrt(X11);
        
        % Scaling variables
        dd = d{i+1};gg = wcgain(i);
        
        % For debugging check tvnorm of scaled closed loop
        if DEBUG
            CL1 = CLoop{i};
            if isa(dd,'tvmat'), dd1 = evalt(dd,CL1.Time); else, dd1 = dd; end
            CLscl = blkdiag(dd1,eye(NY-Nv-Ny))*CL1*blkdiag(1/dd1,eye(NU-Nw-Nu)/gg);
            tvnCLs{i} = tvnorm(CLscl,NE,tvnopt);
            fprintf(' SclPlant CL (TVN): %.4f\n',tvnCLs{i}(2));
        end
        
        %% ScaleOpenLoop for next iteration
        if isa(dd,'tvmat'), dd = evalt(dd,Gunc.Time); end
        Gscl{i+1} = blkdiag(dd,eye(NY-Nv))*Gunc*blkdiag(1/dd,eye(NU-Nw-Nu)/gg,eye(Nu));
        
        %% Display Summary and store info
        if DispFlag
            fprintf(' Summary: Synthesis Gain = %.4f, IQC Analysis Gain = %.4f\n',gs(i),wcgain(i));
        end
        
        % Store info
        if skipinfo
            info{i} = struct('Ks',Ks{i},'CLs',CLs{i},'gs',gs(i),'wcgain',wcgain(i),'d',d{i+1});
        else
            info{i} = struct('Ks',Ks{i},'CLs',CLs{i},'gs',gs(i),'infoK',infoK{i},...
                'wcgain',wcgain(i),'wcinfo',wcinfo(i),'d',d{i+1});
        end
        
        %% Stopping Condition
        if TermDKIterFlag
            if i>2
                % If no significant improvement or increase in gain then abort
                if (wcgain(i-1)-wcgain(i)) <= wcgain(i)*RelTol + AbsTol
                    if wcgain(i-1) < wcgain(i)
                        if Count > 3
                            fprintf(' IQC Analysis Gain has been increasing for last 3 iterations.\n');
                        else
                            Count = Count + 1;
                            fprintf(' IQC Analysis Gain is increasing, Count = %d.\n',Count);
                            continue;
                        end
                    else
                        fprintf(' Stopping tolerance for DK-iteration satisfied.\n');
                    end
                    fprintf(' Terminating iterations...\n');
                    break;
                end
            end
        end
    end
    
catch ME
    fprintf('Unexpected error occurred.\n');
    disp(ME);
    info{i}.ME = ME;
    AbortFlag = true;
end

% Run Finite Number of iterations and select the best controller
% NOTE: Do not count the first iteration
[g,mid] = min(wcgain(1:i));
K = Ks{mid};
CL = CLoop{mid};

% MAX Iter Reached.
if i == NIterMax
    if DispFlag
        fprintf(' Maximum number of iteration reached.\n');
    end
end

% Final Result
if DispFlag && ~AbortFlag
    fprintf(' ---------------------------------------------------------------------------------\n');
    fprintf(' Final Results: WCGain: %.4f, TotalIter: %d, MinIter: %d, TotalTime: %.4f\n',g,i,mid,toc(t0));
    fprintf(' ---------------------------------------------------------------------------------\n');
end
end