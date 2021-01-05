function [A,B,C,D,TS] = ssdata(G)

% XXX - This version does not implement all the functionality
% of the SSDATA command for SS arrays.

GData = G.Data;
[A,B,C,D,TS] = ssdata(GData);
if G.isTimeInvariant
    A = tvmat(A);A.Ts = TS;
    B = tvmat(B);B.Ts = TS;
    C = tvmat(C);C.Ts = TS;
    D = tvmat(D);D.Ts = TS;
else
    Time = G.Time;
    IM = G.InterpolationMethod;
    if isequal(TS,0)
        A = tvmat(A,Time,IM);
        B = tvmat(B,Time,IM);
        C = tvmat(C,Time,IM);
        D = tvmat(D,Time,IM);
    else
        A = tvmat(A,Time,TS);
        B = tvmat(B,Time,TS);
        C = tvmat(C,Time,TS);
        D = tvmat(D,Time,TS);
    end
end
