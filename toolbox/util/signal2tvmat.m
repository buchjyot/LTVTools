function out = signal2tvmat(in)
%% Returns TVMAT for Simulink.SimulationData.Signal input
if isa(in,'Simulink.SimulationData.Signal')
    if length(in.Values.Time) > 2
        out = tvmat(in.Values.Data,in.Values.Time);
    else
        out = tvmat;
    end
else
    error('Input must be of type Simulink.SimulationData.Signal.');
end
end