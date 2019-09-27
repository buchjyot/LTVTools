classdef tltvlib_tvmat < matlab.unittest.TestCase
    
    properties
        testDir = pwd;
        tempFolder;
        model = 'mltvlib_tvmat.slx';
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
        function test1Dtvmat(testCase)
            % Time-Varying gain matrix
            evalin('base','t = 0:0.01:10;');
            evalin('base','K = tvmat(sin(t),t);');
            
            % load the model
            load_system(testCase.model);
            
            % Simulate the model
            sim(testCase.model,0:0.01:10);
            
            % Verify Results
            testCase.verifyEqual(simout,sin(0:0.01:10)','AbsTol',1e-8,...
                'Time-Varying Gain results do not match');
            
            % Clear base workspace
            evalin('base','clear t K');
        end
        
        function testGeneraltvmat(testCase)
            % Time-Varying gain matrix
            evalin('base','t = 0:0.01:10;');
            evalin('base','K = [tvmat(sin(t),t) 1 tvmat(3*t,t); 2 tvmat(cos(t),t) 4];');
            
            % load the model
            load_system(testCase.model);
            
            % Simulate the model
            sim(testCase.model,0:0.01:10);
            simout = reshapedata(tvmat(simout,tout)); %#ok<NODEF>
            t = tout;
            
            % Verify Results
            ONEs = ones(length(t),1);
            expOut = [sin(t) ONEs 3*t 2*ONEs cos(t) 4*ONEs];
            testCase.verifyEqual(simout,expOut,'AbsTol',1e-8,...
                'Time-Varying Gain results do not match');
            
            % Clear base workspace
            evalin('base','clear t K');
        end
    end
end