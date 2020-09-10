function [glb,dwc,info] = tvpowerit(G,Ga,d,CTf,T0,Tf,NL2,x0,MaxIter,StopTol,Display,odeOpt)
%% Power iteration lower bound
% This is an engine function to do power iterations on TVSS objects
U      = cell(MaxIter,1);
Y      = cell(MaxIter,1);
Perf   = zeros(MaxIter,1);
tComp  = zeros(MaxIter,1);
U{1}   = d;
t0     = tic;

% Main for loop
for i = 1:MaxIter
    
    % Start timing
    t1 = tic;
    
    % State Equation
    Y{i} = tvlsim(G,U{i},[T0,Tf],x0,odeOpt);
    
    % Evaluate Performance
    thisY    = Y{i};
    yL2      = thisY(1:NL2,1);
    yE       = tvsubs(thisY(NL2+1:end,1),Tf);
    nyL2     = tvnorm(yL2);
    nyE      = norm(yE);
    Perf(i)  = sqrt(nyE^2 + nyL2^2);
    
    % Display
    if Display
        fprintf(' Iter: %d,\t Perf: %.3f\n',i,Perf(i));
    end
    
    % Stopping Condition
    if (i > 1) && (i < MaxIter)
        % Check if performance that we care is within tolerance
        if Perf(i)-Perf(i-1) <= StopTol*max(Perf([i i-1]))
            break;
        end
    end
    
    % Simulate Costate Equation
    U{i+1} = tvlsim(Ga,yL2,[Tf,T0],CTf'*yE,odeOpt);
    
    % Alignment Condition for U
    U{i+1} = U{i+1}/tvnorm(U{i+1});
    
    % Record single iteration time
    tComp(i) = toc(t1);
end

% Final outputs
tTotal    = toc(t0);
[glb,mid] = max(Perf(1:i));
dwc       = U{mid};
info      = struct('TotalTime',tTotal,'TotalIter',i,'Perf',Perf(1:i),'U',U(1:i),'Y',Y(1:i),'AllIterTime',tComp(1:i));
end