function out = signal2ts(in)
%% Helper function to obtain timeseries object from Simulink.SimulationData.Signal 
if isa(in,'Simulink.SimulationData.Signal')
    out = in.Values;  
else
    error('Input must be of type Simulink.SimulationData.Signal.');
end
end