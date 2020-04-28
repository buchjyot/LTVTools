function n = tvnorm(V,p)
% tvnorm   Lp norm of time-varying signal
%
% For an N-by-1 or 1-by-N TVMAT V, tvnorm(V,p) returns the Lp norm:
%     ( int_0^T (||V(t)||_p)^p dt )^(1/p)
% p can be a positive integer  or inf.
%
% tvnorm(V) is the same as tvnorm(V,2);

% XXX Handle matrices?

%% Input Processing
szV = size(V);
if ~ismatrix(V) || all( szV(1:2) > 1 )
    error('Input to TVNORM must be an N-by-1 or 1-by-N vector')
end
if nargin==1
    p = 2;
end
if isequal(p,-inf)
    error('p must be a positive integer or inf.')
end

%% Compute integrand and integrate over time
I = norm(V,p);
if isinf(p)
    I = max(I);
    n = max(I.Data(:));
else
    if V.isTimeInvariant
        if I.Data==0
            n = 0;
        else
            n = inf;
        end
    else
        I =  I^p;
        n = trapz(I.Time, I.Data(:))^(1/p);
    end
end
