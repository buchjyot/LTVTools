classdef STANDARD_TEST_TEMPLATE < matlab.unittest.TestCase
    
    properties
    end
    
    methods(TestClassSetup)
        function classSetup(testCase)
        end
    end
    
    methods(TestClassTeardown)
        function classCleanup(testCase) %#ok<*MANU>
        end
    end
    
    methods(TestMethodSetup)
        function methodSetup(testCase)
        end
    end
    
    methods(TestMethodTeardown)
        function methodCleanup(testCase)
        end
    end
    
    methods(Test)
        function testMethod1(testCase)
            
        end
        
        function testMethod2(testCase)
            
        end
    end
end