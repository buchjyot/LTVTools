load('twoLinkRobot_LQR.mat','Tlqr');
tvOpt = tvnormOptions('Display','off','AbsTol', 5e-3);

%% Bisection Method
[g1,d1,info1] = tvnormBisect(Tlqr,0,tvOpt)

%% Combined Method
[g2,d2,info2] = tvnormCombined(Tlqr,0,tvOpt)