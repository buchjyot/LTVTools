classdef (CaseInsensitiveProperties = true, ...
        TruncatedProperties = true) poweritWatcher
    
    % WATCHWINDOW  Used by WORSTCASE to display current progress.
    %
    % See also WORSTCASE
    
    properties
        Figure
        Axes
        Title
        Line
    end
    
    methods
        
        function W = poweritWatcher(titleString)
            %% Constructor
            % Create a figure
            W.Figure = figure(...
                'WindowStyle','Docked',...
                'KeyPressFcn',{@LOCAL_keypressfcn},...
                'Interruptible','off');
            
            % Store terminate variable to allow user to quit
            % (Changed by keypressfcn)
            setappdata(W.Figure,'terminate',false);
            
            % Create a set of axes
            W.Axes = axes('NextPlot','Add');
            xlabel('Time (sec)','FontSize',14)
            ylabel('Input','FontSize',14)
            box on; grid on;
            W.Axes.GridLineStyle = '--';
            
            % Create a plot title
            if nargin > 0
                W.Title = title(titleString);
            else
                W.Title = title('');
            end
            
            % This field stores the most recently plotted line
            W.Line = [];
        end
        
        function W = closewatch(W)
            %% CLOSEWATCH
            set(W.Figure,'KeyPressFcn','');
        end
        
        function reply = hasterminated(W)
            %% HASTERMINATED
            % Get the terminate signal from figure's appdata
            if ishandle(W.Figure)
                reply = getappdata(W.Figure,'terminate');
            end
        end
        
        function W = updatewatch(W,U,titleString)
            %% UPDATEWATCH
            % Define constants (RGB colors)
            colorGreen  = [0 1 0];
            colorBlue   = [0 0 1];
            
            % Make sure this figure is active
            figure(W.Figure);
            
            % Make sure the X-Limits are correct
            [u,t] = reshapedata(U);
            if isempty(W.line)
                set(W.Axes,'XLim',[t(1),t(end)]);
            else
                % Make the previous line green
                set(W.Line,'color',colorGreen);
            end
            
            % Plot the current input in Blue
            W.Line = plot(t,u,'color',colorBlue,'LineWidth',1.5);
            
            % Display the objective in the title
            set(W.Title,'String',titleString,'FontSize',14);
            drawnow;
        end
    end
end

function LOCAL_keypressfcn(myFigure,event)
%% LOCAL_keypressfcn
if strncmpi(event.Character,'q',1)  && ishandle(myFigure)
    setappdata(myFigure,'terminate',true);
    disp('Termination requested.  Finishing iteration and closing model...')
    drawnow;
end
end