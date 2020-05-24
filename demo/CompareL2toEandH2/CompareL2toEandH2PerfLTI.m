%% CompareL2toEandH2
% This file runs a for loop to compute the L2toE and H2 norm for different
% horizons. At the end, it creates a MAT file for the data, which then can
% be plotted.

%% Which example do you want to run?
plant_case = 1;
switch plant_case
    case 1
        % LTI SISO Example
        msdfcn = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
        G = msdfcn(1,0.2,1);
        
    case 2
        % LTI SIMO Example
        msdfcn = @(m,b,k) ss([0 1;-k/m -b/m],[0;1/m],[1 0],0);
        G = msdfcn(1,0.2,1);
        G.C = eye(2);
        
    case 3
        % LTI MIMO Example
        Nx = 2; Nu = 2; Ny = 4;
        G = rss(Nx,Ny,Nu);
        G.D = 0;
        
    otherwise
        error('Invalid Option');
end

[Ny,Nu] = size(G);
Nx = order(G);
NE = Ny;
r = rank(G.C);

%% Options
% Horizon
T0 = 0;
Ts = 1e-2;
Tf = [0.01,0.1,0.5,1,2,5,8,10,20,30,50,70];
nT = length(Tf);

% Options
tvnopt = tvnormOptions('OdeSolver','ode45','RelTol',1e-3,'AbsTol',1e-4,'Bounds',[0 10]);
tvopt = tvodeOptions('OdeSolver','ode45');
tvspt = tvlsimOptions('OdeSolver','ode45');

%% Main For Loop
% Memory Allocation
dWc = cell(nT,1);
tvn = cell(nT,1);
Pi = cell(nT,1);
gH2 = zeros(nT,1);
t1 = zeros(nT,1);
t2 = zeros(nT,1);
gE = cell(nT,1);
H2Info = cell(nT,1);
L2Einfo = cell(nT,1);

% For loop for computation
parfor i = 1:nT
    % Time-Varying SS
    Pi{i} = evalt(tvss(G),T0:Ts:Tf(i));
    
    % Compute Norms
    tic;[gE{i},dWc{i},L2Einfo{i}] = tvnorm(Pi{i},NE,tvnopt);t1(i) = toc;
    tic;[gH2(i),H2Info{i}] = tvh2norm(Pi{i},NE,tvopt);t2(i) = toc;
    
    % Display
    fprintf('====================================\n');
    fprintf(' Horizon Tf     = %4.3f s\n',Tf(i));
    fprintf(' L2toENormUB    = %4.3f, Computation Time = %4.3f s\n',gE{i}(2),t1(i));
    fprintf(' TerminalH2Norm = %4.3f, Computation Time = %4.3f s\n',gH2(i),t2(i));
    fprintf(' Check Inequality L2toELB <=    H2 <= sqrt(r)*L2toEUB\n');
    fprintf('                    %.3f <= %.3f <= %.3f\n',gE{i}(1),gH2(i),gE{i}(2)*sqrt(r));
end
fprintf('====================================\n');

%% Monte-Carlo Sims
% This section runs monte-carlo sims for specific example on a specific
% horizons and compares the L2toE and H2 norm interpretations.
switch plant_case
    case {2}
        % Normalize worst-case disturbance & simulate
        normU = 2;
        id = 4; % Choosing the horizon of 1 seconds
        U = normU*dWc{id}/tvnorm(dWc{id});
        Ptv = Pi{id};
        [T0MC,TfMC] = getHorizon(Ptv);
        [Y,X] = tvlsim(Ptv,U,tvspt);
        gLinE = norm(tvsubs(Y,Y.Time(end)))/tvnorm(U);
        
        % Unit Variance White Noise Reponse
        nMC = 100;
        Us = cell(nMC,1);
        UsL2 = zeros(nMC,1);
        Ys = cell(nMC,1);
        Xs = cell(nMC,1);
        fprintf('Running %d Monte-Carlo Sims...\n',nMC);
        parfor j = 1:nMC
            Us{j} = tvmat(randn(TfMC*(1/Ts)+1,1)/sqrt(Ts),T0MC:Ts:TfMC);
            [Ys{j},Xs{j}] = tvlsim(Ptv,Us{j},tvsopt);
        end
        fprintf('Completed %d Monte-Carlo Sims.\n',nMC);
end

%% Save Data
save([mfilename sprintf('_Ex%d',plant_case)]);

%% Plot Data
switch plant_case
    case {1,3}
        % Plot Gain vs Horizons
        figure;clf;
        plot(Tf,cellfun(@(x) x(2), gE),'*b-',Tf,gH2,'rs--');
        xlabel('Horizon (T) (sec)','FontSize',14)
        ylabel('Performance Metric','FontSize',14);
        grid on;box on;
        legend('$\|G\|_{E,[0,T]}$','$\|G\|_{H_2,[0,T]}$','interpreter','latex','Location','southeast','FontSize',14);
        ylim([0 1.8])
        
        % Plot Computational Time
        figure;clf;
        plot(Tf,t1,'*b-',Tf,t2,'rs--');
        xlabel('Horizon (T) (sec)')
        ylabel('Computational Time (sec)');
        legend('$t_{E}$','$t_{H_2}$','interpreter','latex','Location','northwest','FontSize',14);
        grid on;box on;
        
    case 2
        % Gain vs Horizons
        figure;clf;
        plot(Tf,cellfun(@(x) x(2), gE),'*b-',Tf,gH2,'rs--');
        xlabel('Horizon (T) (sec)')
        ylabel('Performance Metric');
        grid on;box on;
        legend('$\|G\|_{E,[0,T]}$','$\|G\|_{H_2,[0,T]}$','interpreter','latex','Location','southeast','FontSize',14);
        ylim([0 2.5])
        
        % Computational Time
        figure;clf;
        plot(Tf,t1,'*b-',Tf,t2,'rs--');
        xlabel('Horizon (T) (sec)')
        ylabel('Computational Time (s)');
        grid on;box on;
        legend('$t_{E}$','$t_{H_2}$','interpreter','latex','Location','northwest','FontSize',14);
        
        % Plot Gaussian ellipsoids as the horizon changes
        X0 = tvsubs(X,0);
        centers = X0';
        legendArray = {};
        figure;clf;hold on;k = 1;
        for i = [3 4 5 6 7]
            plot_gaussian_ellipsoid(centers,H2Info{i}.TerminalVariance,3);
            legendArray{k} = sprintf('T = %.1f',Tf(i)); %#ok<SAGROW>
            k = k + 1;
        end
        % plot(centers(1),centers(2),'ko','MarkerFaceColor','k','MarkerEdgeColor','k');
        grid on;box on;
        xlabel('x_1','FontSize',14);
        ylabel('x_2','FontSize',14);
        legend(legendArray,'location','bestoutside');
        axis equal;
        xlim([-5 5]);ylim([-5 5]);
        
        % Keyboard
        keyboard
        
        % Plot Disks
        figure;clf;hold on;box on; grid on;
        
        % Plot Euclidean Bound
        radii = normU*gE{id}(1);
        angle = linspace(0,2*pi,80);
        p1 = patch(centers(1)+radii*cos(angle),centers(2)+radii*sin(angle),'r');
        p1.EdgeColor = 'r';
        p1.LineWidth = 2;
        alpha(0.1);
        
        % Plot Gaussian Ellipsoid
        plt = plot_gaussian_ellipsoid(centers,H2Info{id}.TerminalVariance,3);
        plt.Color = 'm';
        plt.LineWidth = 2;
        
        % Plot Monte-Carlo Results
        for j = 1:nMC
            Xs0 = tvsubs(Xs{j},0);
            XsTf = tvsubs(Xs{j},Xs{j}.Time(end));
            plot(XsTf(1),XsTf(2),'ko','MarkerFaceColor','g');
        end
        
        % Plot L2toE linearized worst-case disturbace response
        XTf = tvsubs(X,X.Time(end));
        tvplot(X(1),X(2),'b',XTf(1),XTf(2),'bo',...
            'MarkerFaceColor','b','LineWidth',2,'MarkerSize',6);
        
        % Legend
        xlabel('x_1','FontSize',14);
        ylabel('x_2','FontSize',14);
        axis equal;
        
    otherwise
        error('Invalid Option.');
end