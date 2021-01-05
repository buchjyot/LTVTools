function [Ae,Be,Ce,De,CeT,Nz,NeL2,Nd] = extsystem(G,Psi,Nv,Nw,NE)
% This function creates the extended system from G and the IQC multiplier.
% This function assumes the uncertain
%
% Inputs
%   G   - Plant G with uncertainty channels
%   Psi - IQC filter
%   Nv  - Dimentions of Uncertain Output of G
%   Nw  - Dimentions of Uncertain Input to G
%   NE  - Euclidean penalties of interest
% 
% Outputs
%   (Ae,Be,Ce,De)   - State matrices for extended system.
%   (Nz)            - Dimension of IQC filter output, i.e. Nz
%   (NeL2,Nd)       - Dimensions of L2 outputs and disturbance channels

% Input Processing
narginchk(4,5);
if nargin == 4
    NE = 0;
end

%% IQC Multiplier Psi
[Apsi,Bpsi,Cpsi,Dpsi] = ssdata(Psi);
Np = size(Apsi,1);
Bpsi1 = Bpsi(:,1:Nv);
Bpsi2 = Bpsi(:,Nv+1:end);
Dpsi1 = Dpsi(:,1:Nv);
Dpsi2 = Dpsi(:,Nv+1:end);
Nz = size(Psi,1);

%% Extract Euclidean part
NY = size(G,1);
[~,Tf] = getHorizon(G);
idE = NY-NE+1:NY;
[Ag,Bg,Cg,Dg] = ssdata(G);
[CTf,DTf] = tvsubs(Cg(idE,:),Dg(idE,:),Tf);

% Removes Euclidean Part from output matrices
Cg(idE,:) = [];
Dg(idE,:) = [];

% Verify there is no feedthrough from d->e for L2toE
if any(any(DTf))
    error('Feedthrough term d->e must be zero for well-posed L2toE norm.');
end

%% Pull apart nominal system L2 parts
Nx = size(Ag,1);

Bg1 = Bg(1:Nx,1:Nw);
Bg2 = Bg(1:Nx,Nw+1:end);
Cg1 = Cg(1:Nv,1:Nx);
Cg2 = Cg(Nv+1:end,1:Nx);

Dg11 = Dg(1:Nv,1:Nw);
Dg12 = Dg(1:Nv,Nw+1:end);
Dg21 = Dg(Nv+1:end,1:Nw);
Dg22 = Dg(Nv+1:end,Nw+1:end);

[NeL2,Nd] = size(Dg22);

%% Extended System
Ae = [Apsi Bpsi1*Cg1; zeros(Nx,Np) Ag];
Be = [Bpsi1*Dg11+Bpsi2 Bpsi1*Dg12; Bg1 Bg2];
Ce = [Cpsi Dpsi1*Cg1; zeros(NeL2,Np) Cg2];
De = [Dpsi1*Dg11+Dpsi2 Dpsi1*Dg12; Dg21 Dg22];
CeT = [zeros(NE,Np) CTf];
