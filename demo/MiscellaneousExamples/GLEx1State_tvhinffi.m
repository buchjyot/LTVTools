%% HINFFI SISO 1 State Example
% Green, M., & Limebeer, D. J. (2012). Linear robust control. Courier
% Corporation.,
%
% Example 6.2.1

%% Set Options
NE = 0;
RelTol  = 1e-3;
AbsTol  = 1e-5;
Display = 'off';
OdeSolver = 'ode23s';

% Specify tvhinfsynOptions
tvhopts = tvhinfsynOptions('Bounds',[0 1e4],'Display',Display,'RelTol',RelTol,...
    'AbsTol',AbsTol,'OdeSolver',OdeSolver);

% Specify tvnormOptions
tvnopts = tvnormOptions('Bounds',[0 1e4],'Display',Display,'RelTol',RelTol,...
    'AbsTol',AbsTol,'OdeSolver',OdeSolver);

%% System Specifications

% Horizon
T0 = 0;
Tf = [1,2,3,4,5,8,10,20];

% System Matrices & State Space Model
A = 0;
B1 = 1;
B2 = 1;
B = [B1 B2];
C1 = 1;
C = [C1;0];
D11 = [0;0];
D12 = [0;1];
D = [D11 D12];
Gss = ss(A,B,C,D);

% IO Dimentions
Nu = 1;

% Convert LTI to LTV
Gtvss = tvss(Gss);

%% Infinite Horizon FI Controller
if isequal(NE,0)
    [Kih,CLih,gih] = hinffi(Gss,1);
    fprintf(' Inifinte Horizon HINFFI Gain:')
    disp(gih)
    fprintf(' Inifinte Horizon Closed Loop Hinfnorm:');
    disp(hinfnorm(CLih));
    
    % Infinite Horizon L2toL2 Problem converges to gamma = 1 (Analytical Result)
    Tfih = 100;
    Ginih = evalt(Gtvss,T0:0.1:Tfih);
    [tvKih,tvCLih,tvgih,infoih] = tvhinffi(Ginih,Nu,NE,tvhopts);
    fprintf(' Inifinte Horizon Closed Loop L2toL2 Norm:');
    tvnih = tvnorm(tvCLih,NE,tvnopts);
    disp(tvnih);
end

%% Finite Horizon FI Controller
gfh = zeros(length(Tf),1);
for i = 1:length(Tf)
    fprintf(' =============\n');
    fprintf(' Tf = %d\n',Tf(i))
    fprintf(' =============\n');
    
    % Synthesis
    G = evalt(Gtvss,T0:0.01:Tf(i));
    [Kfh,CLfh,gfh(i),infofh] = tvhinffi(G,Nu,NE,tvhopts);
    fprintf(' Synthesis Bisection Gain:')
    disp(gfh(i));
    
    % Verify Bisection Results
    [gtvn,d,info] = tvnorm(CLfh,NE,tvnopts);
    fprintf(' TVNORM Closed Loop Bisection Gain:');
    disp(gtvn);
    
    % Compute L2toE tvnorm
    fprintf(' L2toE Closed Loop:')
    tvnE = tvnorm(CLfh(1,:),1,tvhopts);
    disp(tvnE);
    
    if isequal(NE,0)
        % Plot Riccati Equation solution with Analytical solution
        RiccatiSolP = infofh.Upper.P;
        t = RiccatiSolP.Time;
        
        figure;clf;
        plot(t,RiccatiSolP,'Linewidth',2);hold on;
        T = RiccatiSolP.Time(end);
        phi = sqrt(gfh(i)^-2 - 1); % Analytical solution
        RiccatiSolPAnal = tvmat(phi\tan(-phi*(t-T)),t);
        plot(t,RiccatiSolPAnal,'--r','Linewidth',2);
        legend('Computed Sol.','Analytical Sol.');
        title('Solution to Riccati Diffrential Equation');
        
        % Finite Horizon Analytical Gamma
        gammaAnal = 1/sqrt((pi/(2*T))^2 + 1);
        fprintf(' Finite Horizon Analytical Gamma:')
        disp(gammaAnal);
        
        % Relative Error of Riccati Sol
        figure;clf;
        plot(RiccatiSolP.Time,(RiccatiSolP - RiccatiSolPAnal)/RiccatiSolPAnal,'Linewidth',2);
        title('Relative Error of Riccati Sol');
        
        % Check DownSampling near 0
        figure;clf;
        plot(Kfh.Time,squeeze( Kfh.Data(1,1,:) ),'b',RiccatiSolP.Time,...
            -RiccatiSolP.Data(:),'r--','Linewidth',2);
        title('Check DownSampling near 0');
        xlim([0 0.2]);
    end
end

%% Analysis
% Plot Data
if isequal(NE,0)
    figure;clf;
    plot(Tf,gih*ones(size(Tf)),'r--','LineWidth',3);
    hold on;grid on;
    plot(Tf,gfh,'b-*','LineWidth',3,'MarkerSize',10);
    xlabel('Horizon (T_f)');
    ylabel('L_2 to L_2 Norm');
    title('Closed Loop Induced L_2 Norm vs Time Horizon T_f');
    legend('Infinite Horizon','Finite Horizon');
end