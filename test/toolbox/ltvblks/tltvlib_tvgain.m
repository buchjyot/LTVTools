classdef tltvlib_tvgain < matlab.unittest.TestCase
    
    properties
        testDir = pwd;
        tempFolder;
        model = 'mltvlib_tvgain.slx';
    end
    
    methods(TestMethodSetup)
        function methodSetup(testCase)
            % Create wroking folder for test
            import matlab.unittest.fixtures.TemporaryFolderFixture;
            testCase.tempFolder = testCase.applyFixture(TemporaryFolderFixture);
            
            % Copy the model file to working folder
            copyfile(fullfile(fileparts(which(mfilename)),testCase.model),...
                testCase.tempFolder.Folder,'f');
            cd(testCase.tempFolder.Folder);
        end
    end
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase)
            % Clean up the workspace
            bdclose('all');
            close('all');
            
            % Go back to the testDir and delete the tempFolder by Force
            cd(testCase.testDir);
            [success,msg,msgid] = rmdir(testCase.tempFolder.Folder,'s');
            if ~success
                error(msgid,msg);
            end
        end
    end
    
    methods(Test)
        function testTVGain(testCase)
            % ModelParams
            evalin('base','t = 0:0.01:10;');
            evalin('base','K = [tvmat(cot(t),t) 1; 2 3];');
            evalin('base','inMAT = [tvmat(sin(t),t) 1 0;0 2 3];');
            
            % Clear base workspace
            evalinBase = @(command) evalin('base',command);
            testCase.addTeardown(evalinBase,'clear t K inMAT');
            
            % load the model
            load_system(testCase.model);
            
            % Expected output
            t = 0:0.01:10;
            expout = [tvmat(cos(t),t) tvmat(cot(t),t)+2 3;2*tvmat(sin(t),t) 8 9];
            
            % Simulate the model
            testCase.assumeFail(['Error in ''mltvlib_tvgain/Time-Varying Gain/Time-Varying Matrix/Interpreted',...
                'MATLAB Function''. Evaluation of expression resulted in an invalid output.',...
                'Only finite double vector or matrix outputs are supported ']);
            sim(testCase.model,t);
            
            % Verify Results
            testCase.verifyEqual(simout,expout,'AbsTol',1e-8,...
                'Time Varying Gain Results do not match to expected output.');
        end
    end
end