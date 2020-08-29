function fh = tvsplot(V,varargin)
%% TVSPLOT - Time-varying signal subplots
%
% Input must TVMAT object, followed by LineSpec if any.
%
% If input Y is a (nr x nc) TVMAT then this function creates one figure and
% plots the entry of Y on the grid of subplots vs. time
%
% Example:
% t = 0:0.1:10;
% A = tvmat(sin(t),t);
% B = tvmat(cos(t),t);
% C = [A B; -B -A];
% fh1 = tvsplot(C,'r');

% Create new figure
ltvutil.verifyFH(V);
[nr,nc] = size(V);
TU =  V.TimeUnit;
[VData,VTime] = reshapedata(V);
VName = inputname(1);

% Plot
pid = 1; % subplot id
nout = nargout;
for j = 1:nr
    for k = 1:nc
        subplot(nr,nc,pid);     
        if nout>0
            fh = plot(VTime,VData(:,pid),varargin{:});
        else
            plot(VTime,VData(:,pid),varargin{:});
        end
        pid = pid + 1;
        ylabel(sprintf('%s_{%d%d}',VName,j,k));
        if isequal(j,nr)
            xlabel(sprintf('Time (%s)',TU));
        end
    end
end
end