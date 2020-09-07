function out = simout2struct(simout)
%% Convert simulation input to a structure
nS = length(simout);
out(1:Ns) = struct;
for i = 1:nS
    out(i) = logsout2struct(simout(i).logsout);
end
end