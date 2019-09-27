classdef tAPI_TVMAT < LTVTestBaseClass
    %% TVMAT
    % Basic tests on TVMAT operations
    %
    % <Feature> TVMAT </Feature>
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
            Am = randn(2,0);
            A = tvmat(Am);
            testCase.verifyEqual( isempty(A), isempty(Am) );
            
            Am = randn(2,3);
            A = tvmat(Am);
            testCase.verifyEqual( isempty(A), isempty(Am) );
            
            T = 1:7;
            Am = randn(2,3,4,7);
            A = tvmat(Am,T);
            testCase.verifyEqual( length(A), length(Am(:,:,:,1)) );
            testCase.verifyEqual( ndims(A), ndims(Am(:,:,:,1)) );
            
            testCase.verifyEqual( size(A), size(Am(:,:,:,1)) );
            for dim=1:5
                testCase.verifyEqual( size(A,dim), size(Am(:,:,:,1),dim) );
            end
            [sz1s,sz2s]=size(Am(:,:,:,1));
            [sz1,sz2]=size(A);
            testCase.verifyEqual( [sz1 sz2], [sz1s sz2s] );
            [sz1s,sz2s,sz3s,sz4s,sz5s]=size(Am(:,:,:,1));
            [sz1,sz2,sz3,sz4,sz5]=size(A);
            testCase.verifyEqual( [sz1 sz2 sz3 sz4 sz5], [sz1s sz2s sz3s sz4s sz5s] );
        end
        
        function testUnaryOperationsForConstant(testCase)
            % <Comment>  Test unary operations for DataType="Constant" </Comment>
            rng(0);
            Am = randn(2,3);
            A = tvmat(Am);
            
            B1 = abs(A);
            testCase.verifyEqual( B1.Data, abs(Am) );
            
            B2 = cond(A);
            testCase.verifyEqual( B2.Data, cond(Am) );
            
            B3 = transpose(A);
            testCase.verifyEqual( B3.Data, transpose(Am) );
            
            B4 = tril(A);
            testCase.verifyEqual( B4.Data, tril(Am) );
            
            B5 = mean(A,2);
            testCase.verifyEqual( B5.Data, mean(Am,2) );
        end
        
        function testUnaryOpeartionsForGrid(testCase)
            % <Comment>  Test unary operations for DataType="Grid" </Comment>
            NT = 7;
            T = linspace(0,1,NT);
            Am = crandn(3,3,NT);
            A = tvmat(Am,T);
            
            B1 = sin(A);
            testCase.verifyEqual( B1.Data, sin(Am) );
            
            B2 = det(A);
            for i = 1:NT
                testCase.verifyEqual( B2.Data(:,:,i), det(Am(:,:,i)) );
            end
            
            B3 = ctranspose(A);
            for i = 1:NT
                testCase.verifyEqual( B3.Data(:,:,i), Am(:,:,i)' );
            end
            
            B4 = inv(A);
            for i = 1:NT
                testCase.verifyEqual( B4.Data(:,:,i), inv(Am(:,:,i)) );
            end
            
            TOL = 10*max(size(Am(:,:,1))) * eps(norm(Am(:,:,1)));
            B5 = pinv(A,TOL);
            for i = 1:NT
                testCase.verifyEqual( B5.Data(:,:,i), pinv(Am(:,:,i),TOL) );
            end
        end
        
        function testUnaryOpeartionsOnConstantWithVariableOutputs(testCase)
            
            % <Comment>  Test unary ops on Constant with variable outputs </Comment>
            Am = randn(5);
            A = tvmat(Am);
            
            [Um,Tm] = schur(Am);
            [U,T] = schur(A);
            testCase.verifyEqual( U.Data, Um );
            testCase.verifyEqual( T.Data, Tm );
            
            Tm = schur(Am);
            T = schur(A);
            testCase.verifyEqual( T.Data, Tm );
            
            [Um,Tm] = schur(Am,'complex');
            [U,T] = schur(A,'complex');
            testCase.verifyEqual( U.Data, Um );
            testCase.verifyEqual( T.Data, Tm );
            
            Tm = schur(Am,'complex');
            T = schur(A,'complex');
            testCase.verifyEqual( T.Data, Tm );
        end
        
        function testUnaryOperationsOnGridWithVariableOutputs(testCase)
            % <Comment>  Test unary ops on Grid with variable outputs </Comment>
            NT = 5;
            Am = randn(3,3,NT);
            T = linspace(0,1,NT);
            A = tvmat(Am,T);
            
            [U,T] = schur(A);
            for i=1:NT
                [Um,Tm] = schur(Am(:,:,i));
                testCase.verifyEqual( U.Data(:,:,i), Um );
                testCase.verifyEqual( T.Data(:,:,i), Tm );
            end
            
            T = schur(A);
            for i=1:NT
                Tm = schur(Am(:,:,i));
                testCase.verifyEqual( T.Data(:,:,i), Tm );
            end
            
            [U,T] = schur(A,'complex');
            for i=1:NT
                [Um,Tm] = schur(Am(:,:,i),'complex');
                testCase.verifyEqual( U.Data(:,:,i), Um );
                testCase.verifyEqual( T.Data(:,:,i), Tm );
            end
            
            T = schur(A,'complex');
            for i=1:NT
                Tm = schur(Am(:,:,i),'complex');
                testCase.verifyEqual( T.Data(:,:,i), Tm );
            end
            
        end
        
        function testBinaryOpeartionsForConstantOverConstant(testCase)
            
            % <Comment>  Test binary operations for "Constant"/"Constant" </Comment>
            Am = crandn(3,3);
            A = tvmat( Am );
            Bm = crandn(3,4);
            B = tvmat( Bm );
            C = A*B;
            testCase.verifyEqual( C.Data, Am*Bm );
            
            Am = crandn(3,3);
            A = tvmat( Am );
            Bm = crandn(3,3);
            B = tvmat( Bm );
            C = A+B;
            testCase.verifyEqual( C.Data, Am+Bm );
            
            Am = zeros(0,2);
            A = tvmat(Am);
            Bm = crandn(2,3);
            B = tvmat( Bm);
            C = A*B;
            testCase.verifyEqual( C.Data, Am*Bm );
        end
        
        function testBinaryOpeartionsForConstantOverGrid(testCase)
            % <Comment>  Test binary operations for "Constant"/"Grid" </Comment>
            Am = crandn(3,3);
            A = tvmat( Am );
            NT = 5;
            Bm = crandn(3,3,NT);
            T = linspace(0,1,NT);
            B = tvmat(Bm,T);
            C = A-B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), Am-Bm(:,:,i) );
            end
            
            C = B/A;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), Bm(:,:,i)/Am );
            end
            
            C = B;
            C([1 3],1:2) = 5;
            for i = 1:NT
                testCase.verifyEqual( C.Data([1 3],1:2,i), 5*[1 1; 1 1] );
            end
            
            
            Am = zeros(0,2);
            A = tvmat(Am);
            NT = 5;
            Bm = crandn(2,3,NT);
            T = linspace(0,1,NT);
            B = tvmat(Bm,T);
            C = A*B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), Am*Bm(:,:,i) );
            end
        end
        
        function testBinaryOpeartionsForGridOverGrid(testCase)
            % <Comment>  Test binary operations for "Grid"/"Grid" </Comment>
            NT = 5;
            T = linspace(0,1,NT);
            
            Am = crandn(3,3,NT);
            A = tvmat(Am,T);
            Bm = crandn(3,3,NT);
            B = tvmat(Bm,T);
            
            C = A-B;
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), Am(:,:,i)-Bm(:,:,i) );
            end
            
            C = [A; B; A];
            szC = size(C);
            testCase.verifyEqual( szC, size([Am(:,:,1);Bm(:,:,1); Am(:,:,1)]) );
            for i = 1:NT
                testCase.verifyEqual( C.Data(:,:,i), [Am(:,:,i); Bm(:,:,i); Am(:,:,i)] );
            end
        end
        
        function testRecursiveBinaryOpeartions(testCase)
            % <Comment>  Test recursive binary operations </Comment>
            Am = crandn(3,3);
            A = tvmat( Am );
            Bm = crandn(3,4);
            B = tvmat( Bm );
            C = A*B;
            
            D = [A B C];
            testCase.verifyEqual( D.Data, [Am Bm (Am*Bm)] );
            
            NT = 5;
            T = linspace(0,1,NT);
            Am = crandn(3,3,NT);
            A = tvmat(Am,T);
            Bm = crandn(3,3,NT);
            B = tvmat(Bm,T);
            C = A-B;
            D = blkdiag(A,B,C);
            
            for i = 1:NT
                Ami = Am(:,:,i);
                Bmi = Bm(:,:,i);
                testCase.verifyEqual( D.Data(:,:,i), blkdiag(Ami,Bmi,Ami-Bmi) );
            end
        end
        
        function testLYAPwithConstantDatatype(testCase)
            % <Comment>  Test for LYAP with Constant DataType </Comment>
            Am = randn(4);
            A = tvmat(Am);
            Qm = randn(4);
            Q = tvmat(Qm);
            X = lyap(A,Q);
            testCase.verifyEqual( X.Data, lyap(Am,Qm) );
            
            Bm = randn(4);
            B = tvmat(Bm);
            Cm = randn(4);
            C = tvmat(Cm);
            X = lyap(A,B,C);
            testCase.verifyEqual( X.Data, lyap(Am,Bm,Cm) );
            
            Qm = (Qm+Qm')/2;
            Q = tvmat(Qm);
            Em = randn(4);
            E = tvmat(Em);
            X = lyap(A,Q,[],E);
            testCase.verifyEqual( X.Data, lyap(Am,Qm,[],Em) );
        end
        
        function testLYAPwithGridAndGridOverConstantDatatype(testCase)
            % <Comment>  Test for LYAP with Grid and Grid/Constant DataTypes </Comment>
            NT = 5;
            T = linspace(0,1,NT);
            Am = randn(4,4,NT);
            A = tvmat(Am,T);
            Qm = randn(4,4);
            Q = tvmat(Qm);
            X = lyap(A,Q);
            for i=1:NT
                testCase.verifyEqual( X.Data(:,:,i), lyap(Am(:,:,i),Qm) );
            end
            
            Bm = randn(4,4,NT);
            B = tvmat(Bm,T);
            Cm = randn(4,4,NT);
            C = tvmat(Cm,T);
            X = lyap(A,B,C);
            for i=1:NT
                testCase.verifyEqual( X.Data(:,:,i), lyap(Am(:,:,i),Bm(:,:,i),Cm(:,:,i)) );
            end
            
            Qm = (Qm+Qm')/2;
            Q = tvmat(Qm);
            Em = randn(4,4,NT);
            E = tvmat(Em,T);
            X = lyap(A,Q,[],E);
            for i=1:NT
                testCase.verifyEqual( X.Data(:,:,i), lyap(Am(:,:,i),Qm,[],Em(:,:,i)) );
            end
        end
        
        function testLQRwithConstantDatatype(testCase)
            % <Comment>  Test for LQR with Constant DataType </Comment>
            Am = randn(4);
            A = tvmat(Am);
            Bm = randn(4,2);
            B = tvmat(Bm);
            Qm = randn(4); Qm = Qm*Qm';
            Q = tvmat(Qm);
            Rm = randn(2); Rm = Rm*Rm';
            R = tvmat(Rm);
            
            K = lqr(A,B,Q,R);
            Km = lqr(Am,Bm,Qm,Rm);
            testCase.verifyEqual( K.Data, Km );
            
            [K,S] = lqr(A,B,Q,R);
            [Km,Sm] = lqr(Am,Bm,Qm,Rm);
            testCase.verifyEqual( K.Data, Km );
            testCase.verifyEqual( S.Data, Sm );
            
            [K,S,E] = lqr(A,B,Q,R);
            [Km,Sm,Em] = lqr(Am,Bm,Qm,Rm);
            testCase.verifyEqual( K.Data, Km );
            testCase.verifyEqual( S.Data, Sm );
            testCase.verifyEqual( E.Data, Em );
            
            Nm = randn(4,2);
            N = tvmat(Nm);
            Qm2 = Qm+Nm*inv(Rm)*Nm';  % Ensure Q-N*inv(R)*N>=0
            Q2 = tvmat(Qm2);
            
            K = lqr(A,B,Q2,R,N);
            Km = lqr(Am,Bm,Qm2,Rm,Nm);
            testCase.verifyEqual( K.Data, Km );
            
            [K,S] = lqr(A,B,Q2,R,N);
            [Km,Sm] = lqr(Am,Bm,Qm2,Rm,Nm);
            testCase.verifyEqual( K.Data, Km );
            testCase.verifyEqual( S.Data, Sm );
            
            [K,S,E] = lqr(A,B,Q2,R,N);
            [Km,Sm,Em] = lqr(Am,Bm,Qm2,Rm,Nm);
            testCase.verifyEqual( K.Data, Km );
            testCase.verifyEqual( S.Data, Sm );
            testCase.verifyEqual( E.Data, Em );
        end
        
        function testLQRwithGridAndGridOverConstant(testCase)
            % <Comment>  Test for LQR with Grid and Grid/Constant DataTypes </Comment>
            NT = 4;
            T = linspace(0,6,NT);
            
            Am = randn(4,4,NT);
            A = tvmat(Am,T);
            Bm = randn(4,2,NT);
            B = tvmat(Bm,T);
            Qm = randn(4); Qm = Qm*Qm';
            Q = tvmat(Qm);
            Rm = randn(2); Rm = Rm*Rm';
            R = tvmat(Rm);
            
            K = lqr(A,B,Q,R);
            for i=1:NT
                Km = lqr(Am(:,:,i),Bm(:,:,i),Qm,Rm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
            end
            
            [K,S] = lqr(A,B,Q,R);
            for i=1:NT
                [Km,Sm] = lqr(Am(:,:,i),Bm(:,:,i),Qm,Rm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
                testCase.verifyEqual( S.Data(:,:,i), Sm );
            end
            
            [K,S,E] = lqr(A,B,Q,R);
            for i=1:NT
                [Km,Sm,Em] = lqr(Am(:,:,i),Bm(:,:,i),Qm,Rm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
                testCase.verifyEqual( S.Data(:,:,i), Sm );
                testCase.verifyEqual( E.Data(:,:,i), Em );
            end
            
            
            Nm = randn(4,2);
            N = tvmat(Nm);
            Qm2 = Qm+Nm*inv(Rm)*Nm';  % Ensure Q-N*inv(R)*N>=0
            Q2 = tvmat(Qm2);
            
            K = lqr(A,B,Q2,R,N);
            for i=1:NT
                Km = lqr(Am(:,:,i),Bm(:,:,i),Qm2,Rm,Nm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
            end
            
            [K,S] = lqr(A,B,Q2,R,N);
            for i=1:NT
                [Km,Sm] = lqr(Am(:,:,i),Bm(:,:,i),Qm2,Rm,Nm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
                testCase.verifyEqual( S.Data(:,:,i), Sm );
            end
            
            [K,S,E] = lqr(A,B,Q2,R,N);
            for i=1:NT
                [Km,Sm,Em] = lqr(Am(:,:,i),Bm(:,:,i),Qm2,Rm,Nm);
                testCase.verifyEqual( K.Data(:,:,i), Km );
                testCase.verifyEqual( S.Data(:,:,i), Sm );
                testCase.verifyEqual( E.Data(:,:,i), Em );
            end
        end
        
        function testEVALTwithConstantDatatype(testCase)
            % <Comment>  Test for EVALT with Constant DataType </Comment>
            Am = randn(2);
            A = tvmat(Am);
            NT = 5;
            T = linspace(-2,3,NT);
            B = tvsubs(A,T);
            testCase.verifyEqual( B, repmat(Am,[1 1 NT]) );
            
            B = evalt(A,linspace(-2,3,NT));
            testCase.verifyEqual( B, tvmat(repmat(Am,[1 1 NT]),T) );
        end
        
        function testEVALTwithGridDatatypeLinearInterp(testCase)
            % <Comment>  Test for EVALT with Grid DataType: Linear Interp </Comment>
            AData = cat(3,1,2);
            ATime = [0 5];
            A = tvmat(AData,ATime);
            
            B = tvsubs(A,ATime(1));
            testCase.verifyEqual( B, AData(1,1,1) );
            B = tvsubs(A,ATime(2));
            testCase.verifyEqual( B, AData(1,1,2) );
            B = tvsubs(A,mean(ATime));
            testCase.verifyEqual( B, mean(AData(:)) );
        end
        
        function testEVALTwithGridDatatypeFlatInterp(testCase)
            % <Comment>  Test for EVALT with Grid DataType: Flat Interp </Comment>
            AData = cat(3,1,2);
            ATime = [0 5];
            A = tvmat(AData,ATime,'Flat');
            
            B = tvsubs(A,ATime(1));
            testCase.verifyEqual( B, AData(1,1,1) );
            B = tvsubs(A,ATime(2));
            testCase.verifyEqual( B, AData(1,1,2) );
            B = tvsubs(A,mean(ATime));
            testCase.verifyEqual( B, AData(1,1,1) );
            B = tvsubs(A,0.01*ATime(1)+0.99*ATime(2));
            testCase.verifyEqual( B, AData(1,1,1) );
        end
        
        function testEVALTwithGridDatatypeNearestInterp(testCase)
            % <Comment>  Test for EVALT with Grid DataType: Nearest Interp </Comment>
            AData = cat(3,1,2);
            ATime = [0 5];
            A = tvmat(AData,ATime,'Nearest');
            
            B = tvsubs(A,ATime(1));
            testCase.verifyEqual( B, AData(1,1,1) );
            B = tvsubs(A,ATime(2));
            testCase.verifyEqual( B, AData(1,1,2) );
            B = tvsubs(A,mean(ATime));
            testCase.verifyEqual( B, AData(1,1,2) );
            B = tvsubs(A,0.51*ATime(1)+0.49*ATime(2));
            testCase.verifyEqual( B, AData(1,1,1) );
            B = tvsubs(A,0.49*ATime(1)+0.51*ATime(2));
            testCase.verifyEqual( B, AData(1,1,2) );
        end
        
        function testEVALTwithGridDatatypeSplineInterp(testCase)
            % <Comment>  Test for EVALT with Grid DataType: Spline Interp </Comment>
            T = 0:5;
            NT = numel(T);
            AData = zeros(2,2,NT);
            for i=1:NT
                AData(:,:,i) = [5+7*T(i)^2 -3; -4+T(i)+2*T(i)^3 4*T(i)^3];
            end
            A = tvmat(AData,T,'Spline');
            
            t = linspace(0,5,21);
            Nt = numel(t);
            B = tvsubs(A,t);
            TOL = 1e-8;
            for i=1:Nt
                Ai = [5+7*t(i)^2 -3; -4+t(i)+2*t(i)^3 4*t(i)^3];
                testCase.verifyLessThan( norm( B(:,:,i)- Ai ), TOL );
            end
        end
    end
end