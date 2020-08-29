function r = roots(c)
%% ROOTS without error checking
if c(1)==0
    c = c(find(c~=0,1):end);
end
n = numel(c);
if n>1
    a = diag(ones(1,n-2),-1);
    a(1,:) = -c(2:n) ./ c(1);
    r = eig(a);
else
    r = zeros(0,1);
end
end