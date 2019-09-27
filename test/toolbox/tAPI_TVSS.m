classdef tAPI_TVSS < LTVTestBaseClass
    %% TVSS
    % Basic tests on TVSS operations
    %
    % <Feature> TVSS </Feature>
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
            As = ss(randn(2,0));
            A = tvss(As);
            testCase.verifyEqual( isempty(A), isempty(As) );
            
            As = ss(randn(2,3));
            A = tvss(As);
            testCase.verifyEqual( isempty(A), isempty(As) );
            
            T = 1:7;
            As = rss(2,2,3,4,7);
            A = tvss(As,T);
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
            As = rss(2,2,3,7);
            A = tvss(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,1)) );
            
            T = 1:7;
            As = rss(2,2,3,7);
            A = tvss(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,1)) );
            
            T = 1:7;
            As = rss(2,2,3,4,7);
            A = tvss(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,:,1)) );
            
            T = 1:7;
            As = rss(2,2,3,4,3,7);
            A = tvss(As,T);
            testCase.verifyEqual( ndims(A), ndims(As(:,:,:,:,1)) );
            testCase.verifyEqual( size(A), size(As(:,:,:,:,1)) );
            
        end
        
        function testUnaryOperationsForConstant(testCase)
            % <Comment>  Test unary operations for DataType="Constant" </Comment>
            rng(0);
            As = rss(2,3,4);
            A = tvss(As);
            
            B1 = conj(A);
            testCase.verifyEqual( B1.Data, conj(As) );
            
            B2 = ctranspose(A);
            testCase.verifyEqual( B2.Data, ctranspose(As) );
            
            B3 = norm(A,inf);
            testCase.verifyEqual( B3.Data, norm(As,inf) );
            
            B4 = transpose(A);
            testCase.verifyEqual( B4.Data, transpose(As) );
        end
        
        function testUnaryOperationsForGrid(testCase)
            % <Comment>  Test unary operations for DataType="Grid" </Comment>
            NT = 7;
            T = linspace(0,1,NT);
            As = rss(2,3,3,1,NT);
            A = tvss(As,T);
            
            B1 = conj(A);
            testCase.verifyEqual( B1.Data, conj(As) );
            
            B2 = ctranspose(A);
            for i = 1:NT
                testCase.verifyEqual( B2.Data(:,:,i), ctranspose(As(:,:,i)) );
            end
            
            % State order must be constant across the time grid.
            try
                B3 = inv(A);
                for i = 1:NT
                    testCase.assumeFail( B3.Data(:,:,i), inv(As(:,:,i)) );
                end
            catch ME
                testCase.verifyFalse(isempty(ME),'Exception should not be empty.')
            end
            
            As.D=0;
            A = tvss(As,T);
            B4 = norm(A,2);
            for i = 1:NT
                testCase.verifyEqual( B4.Data(:,:,i), norm(As(:,:,i),2) );
            end
            
            B5 = transpose(A);
            for i = 1:NT
                testCase.verifyEqual( B5.Data(:,:,i), transpose(As(:,:,i)) );
            end
        end
        
        function testBinaryOperationsForConstantOverConstant(testCase)
            % <Comment>  Test binary operations for "Constant"/"Constant" </Comment>
            As = rss(4,3,3);
            A = tvss( As );
            Bs = rss(2,3,4);
            B = tvss( Bs );
            C = A*B;
            testCase.verifyEqual( C.Data, As*Bs );
            
            As = rss(2,3,3);
            A = tvss( As );
            Bs = rss(4,3,3);
            B = tvss( Bs );
            C = A+B;
            testCase.verifyEqual( C.Data, As+Bs );
            
            As = zeros(0,2);
            A = tvss(As);
            Bs = rss(1,2,3);
            B = tvss(Bs);
            C = A*B;
            testCase.verifyEqual( C.Data, As*Bs );
        end
        
        function testBinaryOperationsForConstantOverGrid(testCase)
            % <Comment>  Test binary operations for "Constant"/"Grid" </Comment>
            As = rss(5,3,3);
            A = tvss( As );
            NT = 5;
            Bs = rss(4,3,3,NT);
            T = linspace(0,1,NT);
            B = tvss(Bs,T);
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
                erri = norm(C.Data([1 3],1:2,i)-5*[1 1; 1 1],inf);
                testCase.verifyLessThan( erri, 1e-8 );
            end
            
            As = zeros(0,2);
            A = tvss(As);
            NT = 5;
            Bs = rss(4,2,3,NT);
            T = linspace(0,1,NT);
            B = tvss(Bs,T);
            C = A*B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), As*Bs(:,:,i) );
            end
        end
        
        function testBinaryOperationsForGridOverGrid(testCase)
            % <Comment>  Test binary operations for "Grid"/"Grid" </Comment>
            NT = 5;
            T = linspace(0,1,NT);
            
            As = rss(2,3,3,NT);
            A = tvss(As,T);
            Bs = rss(4,3,3,NT);
            B = tvss(Bs,T);
            
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
        
        function testInterConnectionUtilities(testCase)
            % <Comment>  Test interconnection utilities </Comment>
            Time = linspace(0,5,10)';
            AData = -5+0.1*Time.^2;
            A = tvmat(AData,Time);
            B = 1; C = 1; D=0;
            G = tvss(A,B,C,D);
            M = 2;
            systemnames = 'G M'; %#ok<*NASGU>
            inputvar = '[r]';
            outputvar = '[G]';
            input_to_G = '[r - M]';
            input_to_M = '[G]';
            P = sysic;
            testCase.verifyEqual(P, feedback(G,M) );
            
            C = ss(-1,3,4,0);
            C.u = 'e';
            C.y = 'u';
            G = zpk([],[-1,-1],1);
            G.u = 'u';
            G.y = 'y';
            Sum = sumblk('e = r - y');
            T = connect(G,C,Sum,'r','y');
            Ttv = connect( tvss(G),C,tvss(Sum),'r','y');
            testCase.verifyEqual(Ttv.Data, T );
            
            NT = 4;
            GTime = linspace(0,1,NT);
            GData = rss(3,1,1,NT);
            GData.u = 'u';
            GData.y = 'y';
            
            G = tvss(GData,GTime);
            G.u = 'u';
            G.y = 'y';
            Sum = sumblk('e = r - y');
            Ttv = connect(G,C,Sum,'r','y');
            for i=1:NT
                T = connect( GData(:,:,i),C,Sum,'r','y');
                testCase.verifyEqual(Ttv.Data(:,:,i), T );
            end
        end
    end
end