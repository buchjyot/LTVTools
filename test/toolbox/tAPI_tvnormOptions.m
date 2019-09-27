classdef tAPI_tvnormOptions < matlab.unittest.TestCase
    %% tvnormOptions
    % Basic tests on tvnormOptions and tvodeOptions
    %
    % <Feature> tvnormOptions </Feature>
    % <TestType> API </TestType>
    
    properties(TestParameter)
        ODEsolverToTest = ...
            {'ode45','ode23','ode113','ode23t','ode23s',...
            'ode23tb','ode15s','ode15i'};
    end
    
    methods(Test)
        function testDefaultProperties(testCase)
            % <Comment> Test Default parameters of tvnormOptions </Comment>
            tvnop = tvnormOptions;
            % Update this testpoint if the behavior is changed
            testCase.verifyEqual(numel(properties(tvnop)),6,...
                'There should be 6 properties of tvnormOptions visible to the user')
            testCase.verifyEqual(tvnop.Display,'off',...
                'Default tvnormOptions Display should be "off"');
            testCase.verifyEqual(tvnop.OdeSolver,'ode45',...
                'Default tvnormOptions Solver should be "ode45"');
            testCase.verifyEqual(tvnop.RelTol,1e-3,...
                'Default tvnormOptions RelTol should be "1e-3"');
            testCase.verifyEqual(tvnop.AbsTol,1e-4,...
                'Default tvnormOptions AbsTol should be "1e-4"');
            testCase.verifyEqual(tvnop.AbsTol,1e-4,...
                'Default tvnormOptions AbsTol should be "1e-4"');
            testCase.verifyNotEmpty(tvnop.OdeOptions,...
                'Default tvnormOptions OdeOptions should not be empty');
            OdeOptionsDefault = odeset('RelTol',1e-3,'AbsTol',1e-6);
            testCase.verifyEqual(tvnop.OdeOptions,OdeOptionsDefault,...
                'Default tvnormOptions OdeOptions should be odeset.')
        end
        
        function testOdeSolverSpec(testCase,ODEsolverToTest)
            % <Comment> User can specify OdeSolvers shipping with MATLAB </Comment>
            tvnop = tvnormOptions;
            tvnop.OdeSolver = ODEsolverToTest;
            testCase.verifyEqual(tvnop.OdeSolver,ODEsolverToTest,'Can not set the OdeSolver Property of tvnormOptions');
        end
        
        function testNegOdeSolverSpec(testCase)
            % <Comment> User can not specify any solver other than one shipping with
            % MATLAB </Comment>
            import matlab.unittest.constraints.IssuesWarnings;
            tvnop = tvnormOptions; %#ok<NASGU>            
            testCase.verifyError(@() tvnormOptions('OdeSolver','foo'),...
                'ltvtools:tvodeOptions:OdeSolver');
        end
    end
end