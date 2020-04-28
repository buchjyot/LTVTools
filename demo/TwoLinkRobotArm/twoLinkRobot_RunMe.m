%% References

% [1] Seiler, P., Moore, R.M., Meissen, C., Arcak, M. and Packard, A.,
% 2019. Finite horizon robustness analysis of LTV systems using integral
% quadratic constraints. Automatica, 100, pp.135-143.

% [2] Buch, J., Seiler, P. Finite Horizon Robust Synthesis using IQCs
% submitted to International Journal of Robust and Nonlinear Control, 2020

%% RunMe
% This is a main file that user can run to reproduce the example results
for sequence = 8 % [4,5,6,7,8] % Reproduce results of [2]
    
    switch sequence
        case 0
            % Preliminaries
            twoLinkRobot_BuildLTVModel
            twoLinkRobot_SpecifyOptions
            twoLinkRobot_HinfDesign
            
        case 1
            % Reproduce the results of [1]
            twoLinkRobot_OpenLoop
            twoLinkRobot_LQR
            
        case 2
            % Degradation of Robustness on Finite Horizon
            twoLinkRobot_OpFb
            
            % Plot Results
            twoLinkRobot_OpFb_Plot
            
        case 3
            % Nominal State-Feedback Synthesis
            twoLinkRobot_Hinfsfb
            
        case 4
            % Nominal Output-Feedback Synthesis
            twoLinkRobot_Hinfsyn
            
        case 5
            % Robust Synthesis
            twoLinkRobot_Robsyn
            
        case 6
            % Perform Nominal Analysis on both closed looop
            twoLinkRobot_L2toE_NominalAnalysis
            
            % Perform Robust Analysis on both closed loop
            twoLinkRobot_L2toE_RobustAnalysis
            
        case 7
            % Perform Robust Analysis for different uncertainty level
            twoLinkRobot_UncSweep
            
        case 8
            % Plot results
            twoLinkRobot_NomVsRob_Plot
    end
    
    clearvars -except sequence;
end