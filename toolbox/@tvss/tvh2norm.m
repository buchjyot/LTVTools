function [out,info] = tvh2norm(G,varargin)
%% TVH2NORM Computes finite horizon H2 norm
%
% OUT = TVH2NORM(G) computes finite horizon H2 norm, where G is the
% time-varying state space system.
%
% OUT = TVH2NORM(G,NE) computes combined finite horizon H2 norm over the
% horizon and terminal variance for the NE outputs defined. If NE = NY =
% Total Number of outputs of G, then returned cost is purly terminal
% variance.
%
% [OUT,INFO] = TVH2NORM(G,OPT) computes finite horizon H2 norm of LTV plant
% G on a horizon in which G is defined, OPT is provided as tvodeOptions,
% through which user can specify the ODESOLVER and ODEOPTIONS for CDLE.
%
% [OUT,INFO] = TVH2NORM(G,NE,OPT) computes combined cost, OPT is provided
% as tvodeOptions, through which user can specify the ODESOLVER and
% ODEOPTIONS for CDLE.
%
% Finite Horizon H2 Norm is given by,
%
%           ||G||^2 = (1/T) integral_0toT trace(C(t)*P(t)*C'(t))
%
%           Where, P(t) is a Controllability Gramian satisfying the
%           following lyapunov diffrential equation.
%
%           Pdot(t) = A(t)*P(t) + P(t)*A'(t) + B(t)*B'(t), P(0) = 0
%
% Assumption(s): Feedthrough term D(t) in LTV plant G must be zero.
%                Initial condition variance is assumed to be 0.
%
% Optional input(s): Opt
%
% Example:
% Gss = ss(-1,1,1,0);
% G = tvss(Gss);
% Gt = evalt(G,0:0.01:5);
% norm(Gss,2) % Infinite Horizon 2-Norm
% tvh2norm(Gt) % Finite Horizon 2-Norm
%
% NOTE: Increease the horizon to 50, FH 2-Norm Converges to IH 2-Norm
%
% References:
% Green, M., & Limebeer, D. J. (2012). Linear robust control.
% Courier Corporation. Page 94, Theorem 3.3.1

%% Input Processing
narginchk(1,3);

% Obtain default Opt
[NTf,Opt] = LOCALInputParser(varargin);

% Check Feedthrough condition
[~,~,~,D] = ssdata(G);
if any(D.Data(:))
    warning('ltvtools:tvh2norm:nonzerofeedthrough','The H2 norm is infinite because the system has nonzero feedthrough.');
    out = Inf;
    info = [];
    return;
end

% Check Tspan
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);

% Sizes
NY = size(G,1);
NH2 = 1:NY-NTf;    % Integral cost outputs
NTf = NY-NTf+1:NY; % Terminal cost outputs

% Extract Terminal Part
GTerm = tvss(G.Data(NTf,:),G.Time,G.InterpolationMethod);

% Extract Integral Part
GH2 = tvss(G.Data(NH2,:),G.Time,G.InterpolationMethod);

%% Compute H2 Cost
t1 = tic;
if ~isempty(NTf) && ~isempty(NH2)
    % User wants to compute combined cost
    [IntegralCost,IntegralInfo] = LOCALComputeH2CostOverTheHorizon(GH2,T0,Tf,Opt);
    [TerminalCost,TerminalInfo] = LOCALComputeH2CostAtTerminalTime(GTerm,T0,Tf,Opt);
    % XXX: Consider time normalization for the following equation?
    % i.e. devide by 2?
    out = sqrt((IntegralCost + TerminalCost)/2);
    info = struct('IntegralInfo',IntegralInfo,'TerminalInfo',TerminalInfo);
elseif ~isempty(NTf)
    % User wants to compute Terminal cost
    [TerminalCost,TerminalInfo] = LOCALComputeH2CostAtTerminalTime(GTerm,T0,Tf,Opt);
    out = sqrt(TerminalCost);
    info = TerminalInfo;
elseif ~isempty(NH2)
    % User wants to compute Integral cost
    [IntegralCost,IntegralInfo] = LOCALComputeH2CostOverTheHorizon(GH2,T0,Tf,Opt);
    out = sqrt(IntegralCost);
    info = IntegralInfo;
end
info.TotalTime = toc(t1);
end

function [IntegralCost,IntegralInfo] = LOCALComputeH2CostOverTheHorizon(GH2,T0,Tf,Opt)
%% LOCALComputeH2CostOverTheHorizon
% This function computes the H2 Norm over the horizon and uses an integral
% version of the definition from Green and Limebeer

% Compute Integral of a trace i.e. (1/Tf-T0) integral_0toT trace(C(t)*P(t)*C'(t))
[~,~,CH2,~] = ssdata(GH2);
[GcH2,~,IntegralInfo]  = tvgram(GH2,'c',Opt);
[CH2,GcH2] = evalt(CH2,GcH2,union(CH2.Time,GcH2.Time));
Pt = CH2*GcH2*CH2';
M = trace(Pt);
IntegralCost =  trapz(M.Time,squeeze(M.Data))/(Tf-T0);

% Store info
IntegralInfo.Variance = Pt;  % Output variance as a function of time
IntegralInfo.Gramian = GcH2; % Gramian as a function of time
end

function [TerminalCost,TerminalInfo] = LOCALComputeH2CostAtTerminalTime(GTerm,~,Tf,Opt)
%% LOCALComputeH2CostAtTerminalTime
% This function computes the terminal H2 norm which only uses the final
% time variance.

% Compute Gramian
[~,~,CTerm,~] = ssdata(GTerm);
[GcTerm,~,TerminalInfo]  = tvgram(GTerm,'c',Opt);

% Compute Trace
tUnion = union(CTerm.Time,GcTerm.Time);
[CTerm,GcTerm] = evalt(CTerm,GcTerm,tUnion);
Pt = CTerm*GcTerm*CTerm';
TerminalVariance = tvsubs(Pt,Tf);
TerminalCost = trace(TerminalVariance);

% Store info
TerminalInfo.Variance = Pt;    % Output variance as a function of time
TerminalInfo.Gramian = GcTerm; % Gramian as a function of time
TerminalInfo.TerminalVariance = TerminalVariance;
end

function [NTf,Opt] = LOCALInputParser(V)
%% LOCALInputParser

% Defaults
NTf = 0;
Opt = tvodeOptions;

% Identify distinct input arguments
V1 = V(cellfun(@(x) isa(x,'double'),V));
if ~isempty(V1)
    NTf = V1{1};
end

V2 = V(cellfun(@(x) isa(x,'tvodeOptions'),V));
if ~isempty(V2)
    Opt = V2{1};
end
end