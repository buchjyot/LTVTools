%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2019                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Compute the trajectories for the position and velocity of the Center of
% Mass for a three link planar robot modeling a powered lower limb orthosis
% given the state of the system and its parameters.
% The kinematic equations are published in:
% O.Narvaez-Aroche, A. Packard, and M. Arcak, “Motion planning of the 
% sit-to-stand movement for powered lower limb orthoses”, ASME 2017 Dynamic 
% Systems & Control Conference, October 2017. Tysons, VA, USA.
% http://dx.doi.org/10.1115/DSCC2017-5289
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Function Arguments                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% x: 6 by nt array of System States. 
%   x(1,:): angular position of link 1 relative to x axis (horizontal) [rad].
%   x(2,:): angular position of link 2 relative to link 1 [rad].
%   x(3,:): angular position of link 3 relative to link 2 [rad].
%   x(4,:): angular velocity of link 1 in the inertial frame [rad/s].
%   x(5,:): angular velocity of link 2 [rad/s].
%   x(6,:): angular velocity of link 2 [rad/s].
%
% p: 12 by 1 vector of parameters for the three link robot.
%   p(1): Mass of link 1 [kg]. 
%   p(2): Mass of link 2 [kg]. 
%   p(3): Mass of link 3 [kg]. 
%   p(4): Moment of inertia of link 1 about its CoM [kg.m^2].
%   p(5): Moment of inertia of link 2 about its CoM [kg.m^2].
%   p(6): Moment of inertia of link 3 about its CoM [kg.m^2].
%   p(7): Length of link 1 [m].
%   p(8): Length of link 2 [m].
%   p(9): Length of link 3 [m].
%   p(10): Distance from ankle joint to CoM of link 1 [m]. 
%   p(11): Distance from knee joint to CoM of link 2 [m].
%   p(12): Distance from hip joint to CoM of link 3 [m].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Output                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% y: 4 by nt array of output trajectories.
%   y(1,:): horizontal position of the Center of Mass (CoM) [m].
%   y(2,:): vertical position of the Center of Mass (CoM) [m].
%   y(3,:): horizontal velocity of the Center of Mass (CoM) [m].
%   y(4,:): vertical velocity of the Center of Mass (CoM) [m].
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = x2yThreeLink(x,p)

% Parameters of the system
m1 = p(1);     % Mass of link 1 [kg].
m2 = p(2);     % Mass of link 2 [kg].
m3 = p(3);     % Mass of link 3 [kg].
%I1 = p(4);     % Moment of inertia of link 1 about its CoM.
%I2 = p(5);     % Moment of inertia of link 2 about its CoM.
%I3 = p(6);     % Moment of inertia of link 3 about its CoM.
l1 = p(7);     % Length of link 1 [m].
l2 = p(8);     % Length of link 2 [m].
%l3 = p(9);     % Length of link 3 [m].
lc1 = p(10);   % Distance from ankle joint to CoM of link 1 [m].
lc2 = p(11);   % Distance from knee joint to CoM of link 2 [m].
lc3 = p(12);   % Distance from hip joint to CoM of link 3 [m].

% Constant terms.
k0 = 1/(m1+m2+m3);
k1 = lc1*m1+l1*m2+l1*m3;
k2 = lc2*m2+l2*m3;
k3 = lc3*m3;

% Number of time steps. 
nt = size(x,2);

% Number of outputs. 
ny = 4;

% Array for mapping state to output.
y = zeros(ny,nt);

for i=1:nt
    % Angular positions of the links.
    th = x(1:3,i);
    
    % Sine and cosine functions
    s1 = sin(th(1));
    s12 = sin(th(1)+th(2));
    s123 = sin(th(1)+th(2)+th(3));
    c1 = cos(th(1));
    c12 = cos(th(1)+th(2));
    c123 = cos(th(1)+th(2)+th(3));
    
    % Terms involving sine and cosine functions. 
    k4 = k0*(k1*s1+k2*s12+k3*s123);
    k5 = k0*(k2*s12+k3*s123);
    k6 = k0*k3*s123;
    k7 = k0*(k1*c1+k2*c12+k3*c123);
    k8 = k0*(k2*c12+k3*c123);
    k9 = k0*k3*c123;
    
    % Calculate position coordinates of the CoM.
    y(1,i) = k7;    % x coordinate.
    y(2,i) = k4;    % y coordinate.
    
    % Angular velocities of the links.
    om = x(4:6,i);
    
    % Calculate velocity of the CoM.
    y(3,i) = -om(1)*k4-om(2)*k5-om(3)*k6; % x coordinate.
    y(4,i) = om(1)*k7+om(2)*k8+om(3)*k9;  % y coordinate.
end