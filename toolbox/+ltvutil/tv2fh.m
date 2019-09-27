function Mfh = tv2fh(M)
%% TV2FH Converts TMAT to function_handle 
% Returns function handle for a TVMAT
if isempty(M)
    % Mfh = M;  % Defines Mfh as an empty TVMAT which is slower than
    Mfh = [];
elseif M.isTimeInvariant
    Mfh = @(t) M.Data;
else
    if isequal(M.InterpolationMethod,'Spline') && isempty(M.SplineData)
        M = getSplineData(M);
    end
    Mfh = @(t) tvsubs(M,t);
end