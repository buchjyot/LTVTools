function [out,ndU,ndX0] = stopcond(StopTol,thisU,prevU,thisX0,prevX0)
%% STOPCOND
%
% This function is a general stop condition that is being called in most of
% the solvers used by powerit function.

nin = nargin;
switch nin
    case 3
        tUnion  = union(thisU.Time,prevU.Time);
        [U1,U2] = evalt(thisU,prevU,tUnion);
        Uchange = U1-U2;
        ndU     = tvnorm(Uchange);
        out     = (ndU <= StopTol);
        ndX0    = [];
        
    case 5
        tUnion  = union(thisU.Time,prevU.Time);
        [U1,U2] = evalt(thisU,prevU,tUnion);
        Uchange = U1-U2;
        X0change= thisX0-prevX0;
        ndU     = tvnorm(Uchange);
        ndX0    = norm(X0change);
        out     = (ndU <= StopTol) && (ndX0 <= StopTol);
end
end