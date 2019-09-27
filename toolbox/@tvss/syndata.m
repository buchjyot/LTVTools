function varargout = syndata(G,Ny,Nu,NE)
%% SYNDATA extracts synthesis data from LTV system
%
% G  : TVSS (LTV system)
% Ny : Number of measurements
% Nu : Number of control inputs
% NE : Number of euclidean inputs
%
% Output Order:
% [A,B1,B2,C1,C2,D11,D12,D21,D22,CTf,DTf,Ts,GL2]

% NOTE:
% Ne = NL2 + NE

%% Input Processing
narginchk(3,4);
nin = nargin;
if nin < 4
    NE = 0;
end

%% Extract Eulcidean part if NE~=0
NY= size(G,1);
idE = NY-Ny-NE+1:NY-Ny;
[GL2,CTf,DTf] = LOCALExtractEuclideanData(G,idE);

%% Extract L2 part for H2 and Hinf Synthesis
[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts] = LOCALExtractL2Data(GL2,Ny,Nu);

%% Process outputs
varargout = {A,B1,B2,C1,C2,D11,D12,D21,D22,CTf,DTf,Ts,GL2};
end

%% LOCALExtractEuclideanData
function [G,CTf,DTf] = LOCALExtractEuclideanData(G,idE)
[~,~,CE,DE] = ssdata(tvss(G.Data(idE,:),G.Time,G.InterpolationMethod));
[CTf,DTf] = tvsubs(CE,DE,G.Time(end));

% Removes Euclidean Part from model
G.Data(idE,:) = [];
end

%% LOCALExtractL2Data
function [a,b1,b2,c1,c2,d11,d12,d21,d22,Ts] = LOCALExtractL2Data(G,nY,nU)
% Extract L2 data for H2 or H-infinity synthesis
% Uses Zhou,Doyle,Glover notations
ios = size(G);
nZ = ios(1)-nY;
nW = ios(2)-nU;
idxW = 1:nW;
idxU = nW+1:ios(2);
idxZ = 1:nZ;
idxY = nZ+1:ios(1);
[a,b,c,d,Ts] = ssdata(G);
b1 = b(:,idxW);
b2 = b(:,idxU);
c1 = c(idxZ,:);
c2 = c(idxY,:);
d11 = d(idxZ,idxW);
d12 = d(idxZ,idxU);
d21 = d(idxY,idxW);
d22 = d(idxY,idxU);
end