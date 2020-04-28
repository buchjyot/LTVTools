function [Type,Psi,Psiv,isPsiTV,M] = iqcEngine(Delta)
%% IQC Engine
% This file reads the UserData provided in Delta and returns required
% information needed.
nout = nargout;

% Init
Psi = [];
Psiv = [];
M = [];
isPsiTV = false;

% Extract UserData
[Nw,Nv] = size(Delta);
UserData = Delta.UserData;
if Nv > 1 || Nw > 1
    error('MIMO Delta are not supported yet.');
end

%% Choose a default option for UserData
if isempty(UserData)
    warning('UserData is empty, using the default IQC');
    UserData = [0,-10,1];
end

%% Pack up outputs
switch class(UserData)
    case 'double'
        % Double will be used to provide
        % Type of IQC function that needs to be used
        % Pole and Order of the IQC filter
        Type = UserData(1);
        if nout < 2, return, end
        [Psi,Psiv,isPsiTV] = LOCALIQCFilter(struct('pole',UserData(2),'filterorder',UserData(3)));
        
        % case 'cell'
        % % Cell arrays will be used to specify conic combinations
        % % Provide UserData as {Psi1,M1}, {Psi2,M2} and analysis will solve
        % % for lagrange multiplier lambda
        % Type = 3;
        % if nout < 2, return, end
        % Psi = UserData{1};
        % M = UserData{2};
        % cPsi = class(Psi);
        % if isequal(cPsi,'tvmat') || isequal(cPsi,'tvmat')
        % isPsiTV = true;
        % end
        %
        % case 'struct'
        % % Structures will be used for pre-speciifed IQCs
        % % such as IQC
        % Type = 4;
        % if nout < 2, return, end
        % if UserData{1}.PsiFlag
        % [Psi,Psiv,isPsiTV] = LOCALIQCFilter(UserData{1}.IQCparams);
        % else
        % % Choose Psi to Identity
        % Psi = eye(Nw+Nv);
        % isPsiTV = false;
        % end
        % M = UserData{1}.IQCfunction;
        
    otherwise
        error('Invalid UserData structure specified.')
end
end

function [Psi,Psiv,isPsiTV] = LOCALIQCFilter(params)
%% IQC Multiplier Psi
% Here assumed to be block diagonal
p = params.pole;
q = params.filterorder;

% Dynamics of filter Psi = blkdiag(Psiv,Psiv);
Psi1 = ss(p,1,1,0);
Psiv = ss(1);
for i1=1:q
    Psiv = [eye(i1); zeros(1,i1-1) Psi1]*Psiv;
end
Psi = blkdiag(Psiv,Psiv);

% Filter is not time-varying
isPsiTV = false;
end