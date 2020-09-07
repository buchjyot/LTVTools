%% Specify sequence that you want to run
sequence = 1;
switch sequence
    case 1
        % Reproduce ACC paper results
        TwoStateEx
        TwoStateEx_NominalSynthesis
        TwoStateEx_RobustSynthesis
        TwoStateEx_UncSweep
        TwoStateEx_L2toE_NominalAnalysis
        TwoStateEx_SampleLTIUnc
        TwoStateEx_L2toE_RobustAnalysis
        
    case 2
        % Plot data from MAT file (only)
        TwoStateEx_PlotData
        
    otherwise
        error('Invalid Option');
end