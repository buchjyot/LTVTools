%% Horizon
T0 = 0;
Tf = 3;
Ts = 0.01;
Tgrid = T0:Ts:Tf;

%% System
A = [-0.6 -1;1 0];
B = eye(2);
C = eye(2);
D = 0;
Gss = ss(A,B,C,D);
G = tvss(Gss,Tgrid);
[Ny,Nu] = size(G);

%% Dynamic Response
response_case = 'impulse';

switch response_case
    case 'step'
        % LTI Step Reponse
        Opt1 = stepDataOptions('StepAmplitude',1);
        [y,t,x] = step(Gss,Tf,Opt1);
        
        
        % In LTVTools we also allow user to specify when exactly the step is being
        % applied. This is specified by setting StepTime in the tvsteopOptions
        Opt2 = tvstepOptions('StepAmplitude',1,'StepTime',0);
        [Y,X] = tvstep(G,Tf,Opt2);
        
        % Title String
        titleStr = 'Step Response';
        
    case 'impulse'
        % LTI Impulse Reponse
        [y,t,x] = impulse(Gss,Tf);
        
        % LTV Impulse Response
        Opt = tvodeOptions('OdeSolver','ode23s');
        [Y,X] = tvimpulse(G,Tf,Opt);
        
        % Title String
        titleStr = 'Impulse Response';
        
    case 'initial'
        % Initial Conditions
        x0 = [1;2];
        
        % LTI Initial Condition Reponse
        [y,t,x] = initial(Gss,x0,Tf);
        
        % LTV Initial Condition Reponse
        [Y,X] = tvinitial(G,x0);
        
        % Title String
        titleStr = 'Initial Condition Response';
end

%% Plot
switch response_case
    case {'step','impulse'}
        f1 = figure;
        k = 1;
        for i = 1:Ny
            for j = 1:Nu
                subplot(Ny,Nu,k);
                k = k + 1;
                
                hold on;grid on;box on;
                plot(t,y(:,i,j),'b','LineWidth',2);
                tvplot(Y{i,j},'r--','LineWidth',2);
                
                if isequal(j,1)
                    ylabel(sprintf('To: Out(%d)',i),'FontWeight','normal','Color',[0.5,0.5,0.5]);
                end
                xlabel('');
                if isequal(i,1)
                    title(sprintf('From: In(%d)',j),'FontWeight','normal','Color',[0.5,0.5,0.5]);
                end
                hold off;
            end
        end
        ax=axes(f1,'visible','off');
        ax.XLabel.Visible='on';
        xlabel(ax,'Time (seconds)');
        sgtitle(sprintf('%s',titleStr),'FontWeight','bold');
        
    case {'initial'}
        f1 = figure;hold on;grid on;box on;
        plot(t,y,'b','LineWidth',2);
        tvplot(Y,'r--','LineWidth',2);
        xlabel('Time (seconds)');
        title(sprintf('%s',titleStr),'FontWeight','bold');
end