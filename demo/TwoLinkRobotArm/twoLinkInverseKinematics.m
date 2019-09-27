function [etabar,taubar] = twoLinkInverseKinematics(x2,y2,RobotParam)
% Solve inverse kinematics problem for joint angles, joint velocities,
% and torque inputs corresponding to the motion of link 2 end point.

%% Extract Robot Parameters
beta = RobotParam.beta;
delta = RobotParam.delta;
alpha = RobotParam.alpha;
L1 = RobotParam.L1;
L2 = RobotParam.L2;

%% Solve for joint angles
% The joint angles can be derived from the position of the link 2 end
% point via geometry.  There are two solutions (corresponding to +/- th2
% below).  The specific formulae below are derived in notes entitled
% "Robotics: Forward and Inverse Kinematics" posted by Prof. Damian Gordon.
costh2 = (x2^2+y2^2-L1^2-L2^2)/(2*L1*L2);
th2 = -acos( costh2 );
psi = atan2( L2*sin(th2), L1+L2*costh2 );
phi = atan2(y2,x2);
th1 = phi - psi;

%% Remove Jump Discontinuities
% The formula above yield th2 in [0,pi] and (psi,phi) in [-pi,pi]. Unwrap 
% to remove jump discontinuities in the resulting angles (th1,th2).
th1.Data = unwrap(th1.Data);
th2.Data = unwrap(th2.Data);

%% Use spline to evaluate joint velocities and accels
[th1d,th1dd] = tvdiff(th1);
[th2d,th2dd] = tvdiff(th2);
etabar = [th1; th2; th1d; th2d];

%% Solve for torque inputs
M = [alpha + 2*beta*cos(th2), delta + beta*cos(th2);...
    delta + beta*cos(th2),        delta];
H = [-beta*sin(th2)*th2d, -beta*sin(th2)*(th1d + th2d);...
    beta*sin(th2)*th1d,     0];
taubar =  M*[th1dd;th2dd] + H*[th1d;th2d];