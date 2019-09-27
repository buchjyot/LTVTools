function [K,CL,g,info]= tvh2syn(G,Ny,Nu,Opt)
%% TVH2SYN H2 synthesis on a finite horizon
%
%   [K,CL,GAM] = TVH2SYN(G,NMEAS,NCON) calculates the H2-optimal control
%   law u = K y for the LTV plant G with state-space equations:
%
%       dx =  A(t) x +  B1(t) w +  B2(t) u
%        z = C1(t) x + D11(t) w + D12(t) u
%        y = C2(t) x + D21(t) w + D22(t) u
%
%   NMEAS and NCON specify the numbers of measurements y and controls u (u
%   and y must be the last inputs and outputs of G). The controller K
%   minimizes the H2 norm of CL = LFT(G,K), the closed-loop transfer
%   function from disturbance signals w to error signals z.
%
%   [K,CL,GAM,INFO] = TVH2SYN(G,NMEAS,NCON) calculates the H2
%   controller.
%
%   [K,...] = TVH2SYN(G,NMEAS,NCON,OPT) specifies additional options.
%   Use tvodeOptions to create the option set OPT.
%
% Reference:
%
% Zhou, K., Doyle, J. C., & Glover, K. (1996). Robust and optimal control
% (Vol. 40, p. 146). New Jersey: Prentice hall.
% Page 383, Section 4.5 and Page 387, Section 14.7
% NOTE: LTI Infinite Results are translated to LTV results on Finite Horizon

%% Assumptions
%
% (A1) (A,B2) is stabilizable and (C2,A) is detectable
%
% (A2) D12 has full column rank and D21 has full row rank. i.e. D12'*D12 =
% eye(Nu) and D21*D21' = eye(Ny) have full rank for all times of interest.
%
% (A3) K is causal, linear, time-varying controller satisfying the
% objective ||CL||2norm,[T0,Tf] < g on finite horizon [T0 Tf].

%% Input Processing
narginchk(3,4);
nout = nargout;
nin = nargin;

% Set Defaults
switch nin
    case 3
        Opt = tvodeOptions;
end

% Check Opt
if ~(isa(Opt,'tvodeOptions') || (isa(Opt,'struct') &&...
        isfield(Opt,'OdeOptions') && isfield(Opt,'OdeSolver')))
    error('The options must be of type tvodeOptions or a structure with OdeOptions and OdeSolver as fields.')
end

% LOCALValidateHorizon
ltvutil.verifyFH(G);
[T0,Tf] = getHorizon(G);

%% Call LOCALhtv2synEngine for Design
[K,CL,g,info] = ltvutil.tvh2ric(G,Ny,Nu,T0,Tf,nout,Opt);
end