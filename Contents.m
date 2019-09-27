% LTVTools
% A MATLAB ToolBox for Linear Time-Varying Systems
% Version 1.00 (R2016a)
%
% Grid based time dependent matrices and systems
%   tvmat          - Time-varying matrix.
%   tvss           - Time-varying state-space model.
%   tvumat         - Time-varying uncertain matrix.
%   tvuss          - Time-varying uncertain state-space model.
%
% Manipulation of time dependent models
%   tvsplit        - Extract grid-based LTV model data based on specified horizon.
%   evalt          - Evaluate tvobjects on a user-specified time grid.
%   tvsubs         - Substitutes values for time grid.
%   tvmerge        - Merge two time-varying matrices.
%   tvplot         - Plot time dependent matrix data over time.
%   tvdiff         - Differentiate time-varying matrix.
%
% Model order reduction
%   tvgram         - Compute finite-horizon Gramians for TVSS objects.
%   tvcovar        - Output covariance of LTV system driven by white noise
%
% Robustness and worst-case analysis
%   tvnorm         - Bound on induced-L2 norm for TVSS system.
%   tvh2norm       - Computes H2 norm of TVSS system.
%   tvwcgain       - Worst-case gain of an uncertain LTV system.
%
% Finite-Horizon Controller synthesis
%   tvlqr          - Synthesize a time-varying LQR controller.
%   tvlqg          - Synthesize a time-varying LQG controller.
%   tvkalman       - Synthesize a time-varying Kalman filter.
%   tvhinffi       - Synthesize a full-information Hinf controller.
%   tvh2syn        - Synthesize H2 Optimal Controller.
%   tvhinfsyn      - Synthesize an output-feedback Hinf controller.
%   tvhinfsfb      - Synthesize state-feedback Hinf controller 
%   tvsyn          - A general interface to tvhinfsfb (Ny=0) or tvhinfsyn
%   tvrobsyn       - Synthesize Robust Controller using IQC-synthesis.
%
% Specify Options for LTV analysis and synthesis
%   tvOdeOptions     - Specify OdeSolvers and OdeOptions.
%   tvhinfsynOptions - Create an options object for LTV Hinf-synthesis.
%   tvrobsynOptions   - Create an options object for LTV DK-synthesis.
%   tvnormOptions    - Create options for LTV analysis using tvnorm.
%   tvstepOptions    - Create options for LTV step response.
%
% Time-domain analysis
%   tvlsim       - Simulate time dependent forced response of a TVSS.
%   tvstep       - Simulate time dependent step response of a TVSS.
%   tvinitial    - Simulate initial conditions response of a TVSS.
%
% Simulink
%   ltvlib - Simulink block library for LTV models.
%
% Documentation.
%   Type "doc" and choose Supplemental Software for access to user manual.