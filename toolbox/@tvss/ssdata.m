function [A,B,C,D,TS] = ssdata(G)

% XXX - This version does not implement all the functionality
% of the SSDATA command for SS arrays.

GData = G.Data;
if G.isTimeInvariant
    [A,B,C,D,TS] = ssdata(GData);    
    A = tvmat(A);
    B = tvmat(B);
    C = tvmat(C);
    D = tvmat(D);    
else    
    Time = G.Time;
    IM = G.InterpolationMethod;
    [A,B,C,D,TS] = ssdata(GData);
    A = tvmat(A,Time,IM);
    B = tvmat(B,Time,IM);
    C = tvmat(C,Time,IM);
    D = tvmat(D,Time,IM);
end
