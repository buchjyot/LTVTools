function load_mdlrefs(V)
%% Loads all referenced models
%
% Generally its a good practice to load all the referenced models before
% simulation, that way we can get over the overhead in simulation time that
% is due to loading of referenced models.

allMDLref = find_mdlrefs(V);
for i = 1:numel(allMDLref)
    load_system(allMDLref{i});
end
end