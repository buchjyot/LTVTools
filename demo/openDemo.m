function openDemo(E,varargin)
%% openDemo Opens a Demo for modification and execution.
%
% openDemo(E) opens an example, identified by E, in a new folder. If
% the folder already exists, it opens the existing version.
%
% openExample(E,'workDir', WD) opens an example, identified by E, in a
% folder, identified by WD.

% Check input args
narginchk(1,3);

% Parse and validate inputs
validateChar = @(x) validateattributes(x,{'char'},{'nonempty'});
switch nargin
    case 1
        validateChar(E);
        WD = fullfile(userpath,'LTVToolsDemo',E);
    case 3
        validateChar(E);
        validateChar(varargin{1});
        if ~isequal(varargin{1},'WorkDir')
            error('Working directory must be specified by ''WorkDir'' property.');
        end
        % If dir is specified as char then use as it is otherwise use feval
        % so that 'tempdir' type of arguments can be honored
        if ischar(varargin{2})
            WD = varargin{2};
        else
            WD = feval(varargin{2});
        end
    otherwise
        error('Demo name must be specified with optional working Directory as name-value pair');
end

% If folder exists then just open it otherwise create a folder
if ~isequal(exist(WD,'dir'),7)
    mkdir(WD);
else
    % rmdir(WD);
    % mkdir(WD);
end

% Move the demo files there for execution purpose
copyfile(fullfile(ltvroot,'demo',E,'*'),WD,'f');
cd(WD);
end