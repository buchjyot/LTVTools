classdef tltvlib_tvss < matlab.unittest.TestCase
    
    properties
        testDir = pwd;
        tempFolder;
        model = 'mltvlib_tvss.slx';
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
        function testSISOltvSystem(testCase)
            % System Parameters LTI
            evalin('base','Asys = -1;');
            evalin('base','Bsys = 1;');
            evalin('base','Csys = 1;');
            evalin('base','Dsys = 0;');
            evalin('base','x0sys = 1;');
            
            % System Parameters LTV
            evalin('base','Atvsys = tvmat(Asys);');
            evalin('base','Btvsys = tvmat(Bsys);');
            evalin('base','Ctvsys = tvmat(Csys);');
            evalin('base','Dtvsys = tvmat(Dsys);');
            evalin('base','x0tvsys = x0sys;');
            
            % load the model
            load_system(testCase.model);
            
            % Simulate the model
            sim(testCase.model);
            
            % Verify Results
            testCase.verifyEqual(LTI_RESP,LTV_RESP,'AbsTol',1e-8,...
                'LTI and LTV response do not match');
            
            % Clear base workspace
            evalin('base','clear Asys Bsys Csys Dsys x0sys');
            evalin('base','clear Atvsys Btvsys Ctvsys Dtvsys x0tvsys');
        end
    end
end