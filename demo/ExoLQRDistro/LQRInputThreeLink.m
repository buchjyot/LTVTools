%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez Aroche                                                  %
% Berkeley Center for Control and Identification                          %
% Spring 2019                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the torque at the hips, the torque at the shoulders, the 
% horizontal force at the shoulders, and the vertical force at the 
% shoulders applied by a finite time horizon LQR controller to perform a 
% rest-to-rest maneuver with a powered lower limb orthosis between the 
% sitting position defined by zi and the standing position defined by zf. 
% The nominal parameters of the system used by the controller are specified 
% in the vector pnom.
% The complete control algorithm is published in:
% O. Narvaez-Aroche, A. Packard, and M. Arcak, “Finite time robust 
% control of the Sit-to-Stand movement for powered lower limb orthoses,”
% American Control Conference, June 2018, Milwaukee, WI, USA.
% http://dx.doi.org/10.23919/ACC.2018.8431465
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          Function Arguments                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% T: vector with simulation time steps.
%
% X: state value over time steps.
%   X(:,1): angular position of link 1 relative to x axis (horizontal) [rad].
%   X(:,2): angular position of link 2 relative to link 1 [rad].
%   X(:,3): angular position of link 3 relative to link 2 [rad].
%   X(:,4): angular velocity of link 1 in the inertial frame [rad/s].
%   X(:,5): angular velocity of link 2 [rad/s].
%   X(:,6): angular velocity of link 2 [rad/s].
%
% K: m by n by nt array to interpolate LQR gain at time t.
%
% tgrid: 1 by nt array with values of time for which the m by n LQR gains
%   are known.
%
%
% zi: Coordinates on the z-space for sitting position.
%   zi(1): angular position of link 2 relative to link 1 [rad].
%   zi(2): x coordinate of the Center of Mass (CoM) [m].
%   zi(3): y coordinate of the Center of Mass (CoM) [m].
%
% zf: Coordinates on the z-space for standing position.
%   zf(1): angular position of link 2 relative to link 1 [rad].
%   zf(2): x coordinate of the Center of Mass (CoM) [m].
%   zf(3): y coordinate of the Center of Mass (CoM) [m].
%
% pnom: Nominal values of the parameters of the system, known to controller.
%   pnom(1): Mass of link 1 [kg].
%   pnom(2): Mass of link 2 [kg].
%   pnom(3): Mass of link 3 [kg].
%   pnom(4): Moment of inertia of link 1 about its CoM [kg.m^2].
%   pnom(5): Moment of inertia of link 2 about its CoM [kg.m^2].
%   pnom(6): Moment of inertia of link 3 about its CoM [kg.m^2].
%   pnom(7): Length of link 1 [m].
%   pnom(8): Length of link 2 [m].
%   pnom(9): Length of link 3 [m].
%   pnom(10): Distance from ankle joint to CoM of link 1 [m].
%   pnom(11): Distance from knee joint to CoM of link 2 [m].
%   pnom(12): Distance from hip joint to CoM of link 3 [m].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Output                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% u: array with system Inputs
%   u(:,1): torque applied at the joint of the hips [N.m].
%   u(:,2): torque applied at the joint of the shoulders [N.m].
%   u(:,3): horizontal force applied at the joint of the shoulders [N].
%   u(:,4): vertical force applied at the joint of the shoulders [N].
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = LQRInputThreeLink(T,X,K,tgrid,zi,zf,pnom)

% Acceleration of gravity [m/s^2].
g = 9.81;

% Nominal parameters of the system used by the controller.
m1 = pnom(1);     % Mass of link 1 [kg].
m2 = pnom(2);     % Mass of link 2 [kg].
m3 = pnom(3);     % Mass of link 2 [kg].
I1 = pnom(4);     % Moment of inertia of link 1 about its CoM [kg.m^2].
I2 = pnom(5);     % Moment of inertia of link 2 about its CoM [kg.m^2].
I3 = pnom(6);     % Moment of inertia of link 3 about its CoM [kg.m^2].
l1 = pnom(7);     % Length of link 1 [m].
l2 = pnom(8);     % Length of link 2 [m].
l3 = pnom(9);     % Length of link 3 [m].
lc1 = pnom(10);   % Distance from ankle joint to CoM of link 1 [m].
lc2 = pnom(11);   % Distance from knee joint to CoM of link 2 [m].
lc3 = pnom(12);   % Distance from hip joint to CoM of link 3 [m].

% Constant terms with nominal parameters.
k0 = 1/(m1+m2+m3);
k1 = lc1*m1+l1*m2+l1*m3;
k2 = lc2*m2+l2*m3;
k3 = lc3*m3;

% Array for desired trajectories in the z-space.
zd = zeros(9,1);

% Number of time steps in the simulation.
nt = numel(T);

% Final simulation time.
tf = tgrid(end);

% Array for storing inputs. 
u = zeros(nt,4);

% Loop for retrieving inputs.
for j=1:nt
    
    % Value of time. 
    t = T(j);
    
    % Time polynomials for desired rest-to-rest maneuvers.
    phi = -2*(t/tf)^3 + 3*(t/tf)^2;    % Cubic polynomial of time.
    phip = -6*(t^2/tf^3) + 6*(t/tf^2); % First derivative of cubic polynomial.
    phipp = -12*(t/tf^3) + 6/tf^2;     % Second derivative of cubic polynomial.
    
    % Desired positions on the z-space.
    zd(1:3) = zi + (zf-zi)*phi;
    
    % Desired velocities on the z-space.
    zd(4:6) = (zf-zi)*phip;
    
    % Desired accelerations on the z-space.
    zd(7:9) = (zf-zi)*phipp;
    
    % y coordinate of the CoM of the three link robot in the fully standing
    % position.
    yvert = k0*(k1+k2+k3);
    
    % Compute desired coordinates in the theta space.
    if zd(3) == yvert
        % Velocities and accelerations cannot be computed since matrix V
        % becomes singular. Output array with angles of the robot in the
        % standing position and NaN elements in velocities and acceleration
        % entries.
        thd = [pi/2; 0; 0; NaN; NaN; NaN; NaN; NaN; NaN];
    else
        % Compute angle th1.
        th2 = zd(1);
        alpha = pi+th2;
        if zd(2)>=0
            if zd(3)>=0
                beta = atan(zd(3)/zd(2));
            else
                beta = 2*pi+atan(zd(3)/zd(2));
            end
        else
            beta = pi + atan(zd(3)/zd(2));
        end
        % Squared norm of r1+r2.
        c = k0^2*(k1^2+k2^2+2*k1*k2*cos(th2));
        varphi = asin(max(min(k0*k2*sin(alpha)/sqrt(c),1),-1));
        % Norm of the position of the CoM.
        a = norm(zd(2:3));
        phiang = acos(max(min(((k0*k3)^2-c-a^2)/(-2*sqrt(c)*a),1),-1));
        th1 = beta -phiang +varphi;
        
        % Compute angle th3.
        psi = asin(max(min(sqrt(c)*sin(phiang)/(k0*k3),1),-1));
        th3 = psi -th1 -th2 +beta;
        
        % Compute angular velocities om1 and om3.
        om2 = zd(4);
        s12 = sin(th1+th2);
        s123 = sin(th1+th2+th3);
        c12 = cos(th1+th2);
        c123 = cos(th1+th2+th3);
        f1= zd(5:6)-om2*[-k0*(k2*s12+k3*s123); k0*(k2*c12+k3*c123)];
        V = [-zd(3), -k0*k3*s123; zd(2), k0*k3*c123];
        om13 = V\f1;
        
        % Compute angular accelerations al1 and al3.
        al2 = zd(7);
        q = [-k0*(k2*s12+k3*s123); k0*(k2*c12+k3*c123)]*al2;
        q = q-[zd(2),k0*(k2*c12+k3*c123),k0*k3*c123; zd(3),k0*(k2*s12+k3*s123),k0*k3*s123]*[om13(1)^2;om2^2;om13(2)^2];
        q = q-2*om13(1)*om2*[k0*(k2*c12+k3*c123); k0*(k2*s12+k3*s123)];
        q = q-2*(om13(1)+om2)*om13(2)*[k0*k3*c123; k0*k3*s123];
        al13 = V\(zd(8:9)-q);
        
        % Vector of desired coordinates in the theta space.
        thd = [th1; th2; th3; om13(1); om2; om13(2); al13(1); al2; al13(2)];
    end
    
    % System States
    x = X(j,:)';
    th = x(1:3);
    om = x(4:6);
    
    % Three Link Dynamics equation.
    omsq = [om(1); om(1)+om(2); om(1)+om(2)+om(3)].^2;
    
    % Trigonometric functions
    s2 = sin(th(2));
    s3 = sin(th(3));
    s23 = sin(th(2)+th(3));
    c2 = cos(th(2));
    c3 = cos(th(3));
    c23 = cos(th(2)+th(3));
    sv = sin([th(1); th(1)+th(2); th(1)+th(2)+th(3)]);
    cv = cos([th(1); th(1)+th(2); th(1)+th(2)+th(3)]);
    
    % Mass matrix
    M11 = I1 + I2 + I3 + lc1^2*m1 + m2*(l1^2 + 2*l1*lc2*c2 + lc2^2) + m3*(l1^2 + 2*l1*l2*c2 + 2*l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
    M12 = I2 + I3 + lc2*m2*(l1*c2 + lc2) + m3*(l1*l2*c2 + l1*lc3*c23 + l2^2 + 2*l2*lc3*c3 + lc3^2);
    M13 = I3 + lc3*m3*(l1*c23 + l2*c3 + lc3);
    M22 = I2 + I3 + lc2^2*m2 + m3*(l2^2 + 2*l2*lc3*c3 + lc3^2);
    M23 = I3 + lc3*m3*(l2*c3 + lc3);
    M33 = I3 + lc3^2*m3;
    M = [M11,M12,M13;M12,M22,M23;M13,M23,M33];
    
    % Vector of gravity effects
    fg = -g*[k1,k2,k3;0,k2,k3;0,0,k3]*cv;
    
    % Coriolis effects Matrix
    k4 = l1*(k2*s2+k3*s23);
    k5 = k3*l2*s3;
    C = [k4,-k2*l1*s2+k3*l2*s3,-k3*l1*s23-k3*l2*s3;k4,k5,-k5;l1*k3*s23,k5,0];
    
    % Gravity and Coriolis effects
    f = fg-C*omsq;
    
    % Matrix of generalized force.
    Atau = [[0;0;1],-ones(3,1),-[l1,l2,l3;0,l2,l3;0 0 l3]*sv,[l1,l2,l3;0,l2,l3;0 0 l3]*cv];
    
    % Input constraints for control allocation problem.
    lbu = [-Inf; -Inf; -Inf; 0];      % Lower bounds.
    ubu = [Inf; Inf; Inf; Inf];       % Upper bounds.
    
    % Control allocation weights.
    W = diag([1,1,10,1]);
    
    % Vector for computed torque control approach.
    bCTC = M*thd(7:end)-f;
    
    % Solve constrained linear least-squares problem for nominal input.
    options = optimoptions('lsqlin','Algorithm','active-set','Display','off');
    unom = lsqlin(W,zeros(size(W,1),1),[],[],Atau,bCTC,lbu,ubu,[],options);
    
    % Check if time for interpolation of LQR gain is inside the finite horizon.
    if (t < tgrid(1)) || (t > tf)
        error('The requested time for interpolation of the LQR gain is outside the finite horizon.')
    end
    
    % Number of entries in time array.
    nt = numel(tgrid);
    
    % Find index idx of time array tgrid that contains the minimum value that
    % satisfies t <= tgrid(idx).
    for i = 1:nt
        if  t <= tgrid(i)
            idx = i;
            break
        end
    end
    
    % Interpolate Linear Time-Varying LQR gain at time t.
    if t == tgrid(1)
        Kt = K(:,:,idx);
    else
        % Nearest neighbors of t in time array.
        t0 = tgrid(idx-1);
        t1 = tgrid(idx);
        
        % Find the pair of matrices to be used for interpolation.
        K0 = K(:,:,idx-1);
        K1 = K(:,:,idx);
        
        % Time finite difference.
        dt = t1-t0;
        
        % Rate of change of matrix entries.
        dKdt = (K1-K0)/dt;
        
        % Perform linear interpolation.
        B = K0 - dKdt*t0;
        Kt = dKdt*t + B;
    end
    
    % State of LTV LQR.
    dx = x - thd(1:6);
    
    % Input from TV LQR.
    du = -Kt*dx;
    
    % Control input.
    uctrl = unom + du;
    u(j,:) = uctrl';
end