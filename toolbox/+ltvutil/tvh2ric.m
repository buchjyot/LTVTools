function [K,CL,g,info] = tvh2ric(G,Ny,Nu,T0,Tf,nout,cdreOpt)
%% TVH2SYN Source Code
% RDE based H2 Measurement Feedback Synthesis

% Reference: Zhou, Doyle, Glover, Robust and Optimal Control

%% Define Cost Matrices
S = [];
E = [];
info = [];

% D22 will be required later for controller reconstruction
[~,~,~,~,~,~,D12,D21,D22] = syndata(G,Ny,Nu);

% Assumptions/RankConditions checking
ltvutil.hinfrc(D12,D21,Ny,Nu);

% Normalize D12 and D21 matrices
Gnom = G;
[Gscl,r12,r21] = ltvutil.ioscale(G,Ny,Nu);

% Extract Plant Data
[A,B1,B2,C1,C2,~,D12,D21] = syndata(Gscl,Ny,Nu);
[Nx,Ne,Nd] = syniodim(Gscl,Ny,Nu);
[~,B,C,~] = ssdata(Gscl);

% Costs and State Matrices for X-CDRE
Ah = A - B2*D12'*C1;
Eh = blkdiag(eye(Ne-Nu),zeros(Nu));
Qh = C1'*Eh*C1;
Rh = eye(Nd+Nu);

% Costs and State Matrices for Y-CDRE
Aj = A - B1*D21'*C2;
Ej = blkdiag(eye(Nd-Ny),zeros(Ny));
Qj = B1*Ej*B1';
Rj = eye(Ne+Ny);

% Terminal/Initial Condtion
Fh = zeros(Nx);
Fj = zeros(Nx);

%% Solve CDRE for X and Y
if nout~=4
    % CDRE for X solve backwards in time
    X = cdre(Ah,B,Qh,Rh,S,E,Fh,[Tf T0],cdreOpt);
    % CDRE for Y solve forward in time
    Y = cdre(-Aj,C',-Qj,-Rj,S,E,Fj,[T0 Tf],cdreOpt);
else
    % CDRE for X solve backwards in time
    [X,~,Xdot,Xsol] = cdre(Ah,B,Qh,Rh,S,E,Fh,[Tf T0],cdreOpt);
    % CDRE for Y solve forward in time
    [Y,~,Ydot,Ysol] = cdre(-Aj,C',-Qj,-Rj,S,E,Fj,[T0 Tf],cdreOpt);
    
    % Store Info
    info.Xdot = Xdot;info.Xsol = Xsol;
    info.Ydot = Ydot;info.Ysol = Ysol;
end
info.X = X;
info.Y = Y;

% Use Union of time grids
Tgrid = union(union(Gscl.Time,X.Time),Y.Time);
[X2,Y2,Gscl,Gnom,r21,r12] = evalt(X,Y,Gscl,Gnom,r21,r12,Tgrid);
[A,B1,B2,C1,C2,D11,D12,D21] = syndata(Gscl,Ny,Nu);

%% Design [UNIQUE] H2 optimal measurement feedback controller
F2 = -(D12'*C1 + B2'*X2);
L2 = -(B1*D21' + Y2*C2');

% Case when D11~=0
DK = -D12'*D11*D21';
A2hat = A + B2*F2 + L2*C2;

% Controller Matrices
AK = A2hat-B2*DK*C2;
BK = -(L2-B2*DK);
CK = F2 - DK*C2;
K = tvss(AK,BK,CK,DK);

% Back-substitute scaling factors r12 and r21
K = r12*K*r21;

% Wrap into LFT Interconnection
CL = lft(Gnom,K);

% Compute H2norm of closed loop
g = tvh2norm(CL);

% If D22 feedthrough term exists then based on Page 454: Zhou, K., Doyle,
% J., & Glover, K. (1996). Robust and optimal control.
if any(D22.Data(:))
    K = feedback(K,evalt(D22,K.Time));
end
end