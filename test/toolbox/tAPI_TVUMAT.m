classdef tAPI_TVUMAT < LTVTestBaseClass
    %% TVUMAT
    % Basic tests on TVUMAT operations
    %
    % <Feature> TVUMAT </Feature>
    % <TestType> API </TestType>
    
    methods(TestMethodSetup)
        function methodSetup(testCase) %#ok<MANU>
            rng('default');
        end
    end
    
    methods(Test)
        function testReturnDimensions(testCase)
            % <Comment>  Test functions that return dimensions </Comment>
            rng(0);
            As = umat(randn(2,0));
            A = tvumat(As);
            testCase.verifyEqual( isempty(A), isempty(As) );
            
            As = umat(randn(2,3));
            A = tvumat(As);
            testCase.verifyEqual( isempty(A), isempty(As) );
            
            T = 1:7;
            As = umat(randn(2,3,4,7));
            A = tvumat(As,T);
            testCase.verifyEqual( length(A), length(As(:,:,:,1)) );
            
            testCase.verifyEqual( size(A), size(As(:,:,:,1)) );
            for dim=1:5
                testCase.verifyEqual( size(A,dim), size(As(:,:,:,1),dim) );
            end
            [sz1s,sz2s]=size(As(:,:,:,1));
            [sz1,sz2]=size(A);
            testCase.verifyEqual( [sz1 sz2], [sz1s sz2s] );
            [sz1s,sz2s,sz3s,sz4s,sz5s]=size(As(:,:,:,1));
            [sz1,sz2,sz3,sz4,sz5]=size(A);
            testCase.verifyEqual( [sz1 sz2 sz3 sz4 sz5], [sz1s sz2s sz3s sz4s sz5s] );
            
            T = 1:7;
            As = umat( randn(2,3,7) );
            A = tvumat(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,1)) );
            
            T = 1:7;
            As = umat( randn(2,3,7) );
            A = tvumat(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,1)) );
            
            T = 1:7;
            As = umat( randn(2,3,4,7) );
            A = tvumat(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,:,1)) );
            
            T = 1:7;
            As = umat( randn(2,3,4,3,7) );
            A = tvumat(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,:,:,1)) );
        end
        
        function testUnaryOperationsForConstant(testCase)
            % <Comment>  Test unary operations for DataType="Constant" </Comment>
            rng(0);
            As = umat( randn(2,3,4) );
            A = tvumat(As);
            
            B2 = ctranspose(A);
            testCase.verifyEqual( B2.Data, ctranspose(As) );
            
            B3 = transpose(A);
            testCase.verifyEqual( B3.Data, transpose(As) );
        end
        
        function testUnaryOperationsForGrid(testCase)
            % <Comment>  Test unary operations for DataType="Grid" </Comment>
            NT = 7;
            T = linspace(0,1,NT);
            As = umat( randn(3,3,1,NT) );
            A = tvumat(As,T);
            
            
            B2 = ctranspose(A);
            for i = 1:NT
                testCase.verifyEqual( B2.Data(:,:,i), ctranspose(As(:,:,i)) );
            end
            
            B3 = inv(A);
            for i = 1:NT
                testCase.verifyEqual( B3.Data(:,:,i), inv(As(:,:,i)) );
            end
            
            B5 = transpose(A);
            for i = 1:NT
                testCase.verifyEqual( B5.Data(:,:,i), transpose(As(:,:,i)) );
            end
        end
        
        function testBinaryOperationsForConstantOverConstant(testCase)
            % <Comment>  Test binary operations for "Constant"/"Constant" </Comment>
            As = umat( randn(3,3) );
            A = tvumat( As );
            Bs = umat( randn(3,4) );
            B = tvumat( Bs );
            C = A*B;
            testCase.verifyEqual( C.Data, As*Bs );
            
            As = umat( randn(3,3) );
            A = tvumat( As );
            Bs = umat( randn(3,3) );
            B = tvumat( Bs );
            C = A+B;
            testCase.verifyEqual( C.Data, As+Bs );
            
            As = zeros(0,2);
            A = tvumat(As);
            Bs = umat( randn(2,3) );
            B = tvumat(Bs);
            C = A*B;
            testCase.verifyEqual( C.Data, As*Bs );
        end
        
        function testBinaryOperationsForConstantOverGrid(testCase)
            % <Comment>  Test binary operations for "Constant"/"Grid" </Comment>
            As = umat( randn(3,3) );
            A = tvumat( As );
            NT = 5;
            Bs = umat( randn(3,3,NT) );
            T = linspace(0,1,NT);
            B = tvumat(Bs,T);
            C = A-B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), As-Bs(:,:,i) );
            end
            
            C = B/A;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), Bs(:,:,i)/As );
            end
            
            C = B;
            C([1 3],1:2) = 5;
            for i = 1:NT
                testCase.verifyEqual( C.Data([1 3],1:2,i), umat(5*[1 1; 1 1]) );
            end
            
            As = zeros(0,2);
            A = tvumat(As);
            NT = 5;
            Bs = umat( randn(2,3,NT) );
            T = linspace(0,1,NT);
            B = tvumat(Bs,T);
            C = A*B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), As*Bs(:,:,i) );
            end
        end
        
        function testBinaryOperationsForGridOverGrid(testCase)
            % <Comment>  Test binary operations for "Grid"/"Grid" </Comment>
            NT = 5;
            T = linspace(0,1,NT);
            
            As = umat( randn(3,3,NT) );
            A = tvumat(As,T);
            Bs = umat( randn(3,3,NT) );
            B = tvumat(Bs,T);
            
            C = A-B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), As(:,:,i)-Bs(:,:,i) );
            end
            
            C = [A; B; A];
            szC = size(C);
            testCase.verifyEqual( szC, size([As(:,:,1);Bs(:,:,1); As(:,:,1)]) );
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), [As(:,:,i); Bs(:,:,i); As(:,:,i)] );
            end
        end
    end
end