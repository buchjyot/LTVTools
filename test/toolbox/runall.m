function results = runall()
%% Run All Tests
allTests = dir('t*.m');
results = runtests({allTests.name},'verbosity',4);