function results = runAllTests()
%% Run all tests

% Recursively find all runall.m files under testroot and run them
runFiles = dir(fullfile(ltvroot,'test','**/runall.m'));
results = cell(numel(runFiles),1);
passedFlag = true;

% Run all the folders
for i = 1:numel(runFiles)
    fprintf('=========================================================\n');
    fprintf('### Running Test from the following directory\n%s\n',...
        runFiles(i).folder);
    fprintf('=========================================================\n');
    results{i} = fevalin(fullfile(runFiles(i).folder),'runall');
    
    % If all the tests did not pass then raise the flag
    if ~all([results{i}.Passed])
        passedFlag = false;
    end
end

% Print Status if there are failures or not
if ~passedFlag
    fprintf('### There were failures in the test results.\n');
else
    fprintf('### There were no test failures.\n');
end
end

%% Evaluate function at specific location
function out = fevalin(location, fcn)
oldfolder = cd(location);
out = feval(fcn);
cd(oldfolder);
end