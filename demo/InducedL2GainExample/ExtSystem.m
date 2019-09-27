function [Ae,Be,Ce,De,Nz,Ne,Nd] = ExtSystem(G,v,p)
% This function creates the extended system from G and the
% IQC multiplier. This function assumes the uncertain
%
% Inputs
%   G - Nominal LTV system (as a TVSS). The uncertain system is
%       Fu(Delta,G) where Delta is a SISO, LTI, uncertainty
%       satisyfing the unit norm bound ||Delta|| <= 1.
%   (v,p) - The IQC multiplier is blkdiag(Psiv'*X11*Psiv,-Psiv'*X11*Psiv)
%           where Psiv:=[1; 1/(s-p); ...; 1/(s-p)^v] and X11>=0.
%
% Outputs
%   (Ae,Be,Ce,De) - State matrices for extended system.
%   Nz- Dimension of IQC filter output, i.e. Nz=2*(v+1)
%   (Ne,Nd) - Dimensions of error and disturbance channels

%% Pull apart nominal system
[Ag,Bg,Cg,Dg] = ssdata(G);
Nx = size(Ag,1);
Bg1 = Bg(1:Nx,1);
Bg2 = Bg(1:Nx,2:end);
Cg1 = Cg(1,1:Nx);
Cg2 = Cg(2:end,1:Nx);

Dg11 = Dg(1,1);
Dg12 = Dg(1,2:end);
Dg21 = Dg(2:end,1);
Dg22 = Dg(2:end,2:end);

[Ne,Nd] = size(Dg22);

%% IQC Multiplier Psi
Psi1 = ss(p,1,1,0);
Psiv = ss(1);
for i1=1:v
    Psiv = [eye(i1); zeros(1,i1-1) Psi1]*Psiv;
end
Psi = blkdiag(Psiv,Psiv);

[Apsi,Bpsi,Cpsi,Dpsi] = ssdata(Psi);
Bpsi1 = Bpsi(:,1);
Bpsi2 = Bpsi(:,2);
Dpsi1 = Dpsi(:,1);
Dpsi2 = Dpsi(:,2);
Nz = 2*(v+1);

%% Extended System
Np = 2*v;
Ae = [Apsi Bpsi1*Cg1; zeros(Nx,Np) Ag];
Be = [Bpsi1*Dg11+Bpsi2 Bpsi1*Dg12; Bg1 Bg2];
Ce = [Cpsi Dpsi1*Cg1; zeros(Ne,Np) Cg2];
De = [Dpsi1*Dg11+Dpsi2 Dpsi1*Dg12; Dg21 Dg22];
