%% References

% [1] Seiler, P., Moore, R.M., Meissen, C., Arcak, M. and Packard, A.,
% 2019. Finite horizon robustness analysis of LTV systems using integral
% quadratic constraints. Automatica, 100, pp.135-143.

%% RunMe
% This is a main file that user can run to reproduce the example results
sequence = 1;

switch sequence
    case 1 
        % Reproduce the results of [1]
        twoLinkRobot_BuildLTVModel
        twoLinkRobot_OpenLoop
        twoLinkRobot_LQR
end