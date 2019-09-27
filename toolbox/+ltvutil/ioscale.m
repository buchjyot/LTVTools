function [Pscl,r12,r21] = ioscale(P,nY,nU)
%% IOScale
% Scale E and D by unitary, and U and Y by invertible to change D12 and D21
% into the simplified [0;I], and [0 I] form.  Return scaled plant
% invertible transformations on Y and U.

[ap,b1,b2,c1,c2,d11,d12,d21,d22] = syndata(P,nY,nU);
[nE,nD] = size(d11);

[q12,r12] = qr(d12);  % Should be full-column rank
q12 = q12(:,[(nU+1):end 1:nU])';
r12 = eye(nU)/r12(1:nU,:);

[q21,r21] = qr(d21'); % Should be full-column rank
q21 = q21(:,[(nY+1):end 1:nY]);
r21 = eye(nY)/(r21(1:nY,:)');

%  Scale the matrices to:
%    q12*d12*r12 = [0; I] and r21*d21*q21 = [0 I]
c1 = q12*c1;
c2 = r21*c2;
cp = [c1;c2];
b1 = b1*q21;
b2 = b2*r12;
bp = [b1,b2];
d11 = q12*d11*q21;
d12 = [zeros(nE-nU,nU); eye(nU)];  % d12 = q12*d12*r12;
d21 = [zeros(nY,nD-nY), eye(nY)];  % d21 = r21*d21*q21;

% XXX - The formulas derived in the functions above (in the main and other
% supporting functions) assume d22=0.   This is accomplished with an
% initial loop transformation, which is "undone" in the controller
% reconstruction.
dp = [d11 d12;d21 0*d22];
Pscl = tvss(ap,bp,cp,dp);
end