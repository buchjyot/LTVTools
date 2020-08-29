%% PolySysEx
% This demo needs polynomial systems toolbox and multipoly from
% https://dept.aem.umn.edu/~AerospaceControl/

% If you don't have them installed then same example is implemented in
% simulink and can be run by changing the following switch case option.

%% Choose Example
exampleCase = 'simulink';
switch exampleCase
    
    case 'polysys'
        % Create a model of the system.
        % First, |polynomial| variables are created using the |pvar| command.  Then,
        % these variables are used to define the functions |f| and |g|, which are
        % also |polynomial| variables.
        pvar x1 x2 u
        states = [x1;x2];
        inputs = u;
        f = [ -x1 + x2 - x1*x2^2 ; -x2*x1^2 - x2 + u ];
        g = states;
        
        % Then, a |polysys| object is created from the polynomials |f| and |g|.
        sys = polysys(f,g,states,inputs);
        
        % The |polynomial| objects |states| and |inputs| specify the ordering of the
        % variables. That is, by setting |states(1) = x1|, we specify that |f(1)|
        % is the time derivative of |x1|.
        
    case 'simulink'
        % Note: This model is configured to stop at T = 10 seconds.
        sys = 'NLSys';
        load_system(sys);
end

%% Power Iterations
T0 = 0;
Tf = 10;
pSpec = poweritSignalSpec('NE',2,'InitialInput','ones');
pOpt = poweritOptions('Display','on');
[glb,dwc,info] = powerit(sys,[T0,Tf],pSpec,pOpt);

%% Plot
figure(1);clf;
tvplot(dwc,'LineWidth',2);
box on;grid on;
ylabel('Worst-Case Disturbance');
xlabel('Time (sec)');