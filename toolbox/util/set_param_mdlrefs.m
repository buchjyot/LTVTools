function set_param_mdlrefs(V,varargin)
%% This function does set_param to all the referenced model
%
% Example: 
% set_param_mdlrefs('mTop','AbsTol','1e-3');

nin = nargin;

if isequal(nin,0) || isequal(nin,1)
    return
elseif isequal(mod(nin,2),1)
    % Set Name value pairs
    allMDLref = find_mdlrefs(V);
    for i = 1:numel(allMDLref)
        set_param(allMDLref{i},varargin{:});
    end
else
    error('Input must be in the form of name-value pairs.');
end
end