% System
G = tvss(5,1,1,0);

% Power iterations
pOpt = poweritOptions('Display','on');
Tf = [0.1 1 2 5 10 20 30 35];
Nt = length(Tf);
for i = 1:Nt
    H = [0,Tf(i)];
    Gt = evalt(G,H);
    gPower(i) = powerit(Gt,H,pOpt); %#ok<SAGROW>
    tOpt = tvnormOptions('Bounds',[0 gPower(i)*1.5],'Display','on');
    gbnd = tvnormb(Gt,tOpt);
    gRDE(i) = gbnd(2); %#ok<SAGROW>
end

figure,
semilogy(Tf,gRDE,'-bo',Tf,gPower,'--g^','LineWidth',2)
grid on;box on;xlabel('Horizon T (sec.)');
ylabel('Induced L_2 Gain');
legend('RDE Upper Bound','Power Iteration Lower Bound','Location','southeast')