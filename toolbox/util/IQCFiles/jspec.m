function [Psi,Jpq,X] = jspec(Pi)
% function [Psi,M] = jspec(Pi)
%
% This function constructs a J-spectral factorization for Pi. In other
% words, it factorizes Pi as:
%   Pi = Psi'*M*Psi
% where Psi is stable with a stable inverse and M = blkdiag(I,-I).
%
% Input
%  Pi:  Parahermitian IQC multiplier as a TF/SS object.
%
% Outputs
%  Psi: Auxiliary system as an SS object. Psi is stable and has
%          a stable inverse.
%  M: Constant matrix of the form blkdiag(I,-I)
%
% Note: Some technical conditions on Pi are required to obtain a J-spectral
% factorization.  A sufficient condition is for the upper left/lower right
% blocks of Pi to satisfy Pi_11>0 and Pi_22<0 on the imaginary axis. This
% function assumes that the technical conditions required are satisfied by
% the multiplier Pi.

%% Ensure Minimum Realization
% Construct minimal, state-space realization of IQC multiplier
Pi = ss(Pi);
Pi = minreal(Pi);

%% IQC Factorization
% Factorize Pi in to Psi'*M*Psi where Psi is stable.
% Use method described in Section 7.3 in "A Course in Hinf Control
% Theory," by B. Francis.  This method is also described in the proof
% of the first factorization theorem in the IQC paper.
[GS,GU] = stabsep(Pi); %#ok<ASGLU>
nx = size(GS.A,1);
m = size(Pi,1);

Q = zeros(nx);
R = GS.D;
S = GS.C';

Apsi = GS.A;
Bpsi = GS.B;

%% J-spectral Factorization
% R = U*T*U' and T should be real, diagonal since R=R'
% Factor R = W'*Jpq*W where p and q are the # of positive
% and negative evals of R, respectively.
[U,T] = schur(R);
ev = real( diag(T) );
pidx = find( ev>0 );
qidx = find( ev<0 );
W = diag( sqrt([ev(pidx); -ev(qidx)]) );
W = W*U(:,[pidx;qidx])';

p = length(pidx);
q = length(qidx);

if ~isequal(p+q,m)
    error('Psi should not have any eigenvalues on imaginary axis.');
end

%% Process Outputs
Jpq = blkdiag( eye(p), -eye(q) );
X = care(Apsi,Bpsi,Q,R,S);
Psi = ss(Apsi,Bpsi,Jpq*(W'\(Bpsi'*X+S')),W);