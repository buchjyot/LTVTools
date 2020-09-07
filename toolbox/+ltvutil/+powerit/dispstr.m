function out = dispstr(varargin)
%% DISPSTR

% Input processing
nin = nargin;
switch nin
    case 3
        [id,pf,deltaU] = deal(varargin{:});dx0 = [];
    case 4
        [id,pf,deltaU,dx0] = deal(varargin{:});
end

% Prepare outputs
c1 = isempty(deltaU) || isequal(deltaU,0);
c2 = isempty(dx0) || isequal(dx0,0);
switch num2str([c1,c2])
    case '1  1'
        out = sprintf(' Iter: %d,\t Perf: %.3f',id,pf);
    case '0  1'
        out = sprintf(' Iter: %d,\t Perf: %.3f,\t dU: %.3f\t',id,pf,deltaU);
    case '1  0'
        out = sprintf(' Iter: %d,\t Perf: %.3f,\t dx0: %.3f\t',id,pf,dx0);
    otherwise
        out = sprintf(' Iter: %d,\t Perf: %.3f,\t dU: %.3f,\t dx0: %.3f',id,pf,deltaU,dx0);
end
end