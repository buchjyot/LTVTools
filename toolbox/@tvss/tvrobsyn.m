function [K,CL,g,info] = tvrobsyn(Gunc,Delta,Ny,Nu,varargin)
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

%% Initialize & Input Processing
NE = 0;
narginchk(4,6);
nin = nargin;
% nout = nargout;
switch nin
    case 5
        NE = varargin{1};
        Opts = tvrobsynOptions;
    case 6
        % Default Order
        NE = varargin{1};
        Opts = varargin{2};
end

% If NE is empty [] then make it zero
if isempty(NE)
    NE = 0;
end

%% TVROBSYN Engine
% RDE based Synthesis and RDE/DLMI based Analysis
nout = nargout;

% Set Special Flags
ABORTFlag = false;

% Read Options
NIterMax = Opts.MaxIter;
DispFlag = isequal(Opts.Display,'on');
tvhopt = Opts.SynthesisOptions;
tvwcopt = Opts.AnalysisOptions;
DEBUG = Opts.DebugMode;

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
[Nw,Nv] = size(Delta);
[T0,Tf] = getHorizon(Gunc);
Tspan = [T0,Tf];

% User can specify either function handle or double value for this
Dinit = Opts.InitialDScale;
if isa(Dinit,'function_handle')
    Dinit = Dinit(Nw,Nv);
end

% Gscl stands for scalled open-loop plant
[NY,NU] = size(Gunc);
d{1} = Dinit;
Gscl{1} = blkdiag(Dinit,eye(NY-Nv))*Gunc*blkdiag(1/Dinit,eye(NU-Nw));

% Start timing
t0 = tic;
try
    
    %% DK-Iterations
    for i = 1:NIterMax
        
        % Start Display
        if DispFlag
            fprintf(newline);
            fprintf(' ===> Iteration #%d:\n',i);
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
        % XXX : Assumes Delta is same across the time grid.
        GextKt = evalt(Gext,Ks{i}.Time);
        CLoop{i} = lft(GextKt,Ks{i});
        
        % D-Step
        if DispFlag
            fprintf(' ### IQC Analysis Step: \n');
        end
        [wcgain(i),wcinfo{i}] = tvwcgain(CLoop{i},Delta,NE,tvwcopt);
        if DispFlag
            fprintf(' IQC Analysis Gain: = %.4f \n', wcgain(i));
        end
        
        % D*D is equal to Psiv'*X11*Psiv
        X11 = wcinfo{i}.LMI.X11;
        d{i+1} = LOCALSpectralFact(Delta,X11,Tspan);
        
        % Scaling variables
        dd = d{i+1};gg = wcgain(i);sg = sqrt(gg);
        
        % For debugging check tvnorm of scaled closed loop
        % Numerically it is better to have the factor of sqrt(gamma) on the
        % performance channel both inputs and outputs.
        if DEBUG
            CL1 = CLoop{i};
            dd1 = LOCALRegridDScale(dd,CL1.Time);            
            CLscl = blkdiag(dd1,eye(NY-Nv-Ny)/sg)*CL1*blkdiag(1/dd1,eye(NU-Nw-Nu)/sg);
            tvnCLs{i} = tvnorm(CLscl,NE,tvnopt);
            fprintf(' SclPlant CL (TVN): %.4f\n',tvnCLs{i}(2));
        end
        
        %% ScaleOpenLoop for next iteration
        dd = LOCALRegridDScale(dd,Gunc.Time);
        Gscl{i+1} = blkdiag(dd,eye(NY-Nv)/sg)*Gunc*blkdiag(inv(dd),eye(NU-Nw-Nu)/sg,eye(Nu));
        
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
    fprintf(newline);
    fprintf(' ========================================\n');
    fprintf(' Unexpected error occurred.\n');
    for i = 1:length(ME.stack)
        fprintf(' ==> Error in line %d in function %s of the file %s\n',ME.stack(i).line,ME.stack(i).name,ME.stack(i).file);
    end
    fprintf(' ========================================\n');
    throw(ME);
    ABORTFlag = true;
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
if DispFlag && ~ABORTFlag
    fprintf(' ---------------------------------------------------------------------------------\n');
    fprintf(' Final Results: WCGain: %.4f, TotalIter: %d, MinIter: %d, TotalTime: %.4f\n',g,i,mid,toc(t0));
    fprintf(' ---------------------------------------------------------------------------------\n');
end
end

function out = LOCALRegridDScale(in,t)
%% LOCALRegridDScale
switch class(in)
    case {'tvmat','tvss'}
        out = evalt(in,t);
    otherwise
        out = in;
end
end

function out = LOCALSpectralFact(Delta,X11,Tspan)
%% LOCALSpectralFact
FH_SPECTRAL_FACT = true;
[m,n] = size(X11);
DYNAMIC_IQC = (prod(m,n) > 1);

if DYNAMIC_IQC
    % Spectral Factorization
    [~,~,Psiv] = ltvutil.iqcEngine(Delta);
    if FH_SPECTRAL_FACT
        % Finite Horizon
        out = tvspectralfact(Psiv,X11,Tspan);
    else
        % Infinite Horizon
        [Dv,S] = spectralfact(Psiv,X11);
        out = chol(S)*Dv;
    end
else
    out = sqrt(X11);
end
end