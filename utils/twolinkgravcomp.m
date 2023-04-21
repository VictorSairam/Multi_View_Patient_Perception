% Two-link Gravity Compensation Example    
clear all
close all
clc

L1 = 0.3; % [m]
L2 = 0.3; % [m]
m1 = 1;   % [kg]
m2 = 1;   % [kg]

robot = SerialLink([Revolute('d', 0, 'a', L1, 'alpha', 0, 'offset', 0), ...
                    Revolute('d', 0, 'a', L2, 'alpha', 0, 'offset', 0)], ...
                    'name', 'Two-Link Robot');
                
robot.teach(zeros(1,2));

% M matrices
M01 = [1 0 0 L1;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
  
M2 = [1 0 0 L1+L2;
      0 1 0   0;
      0 0 1   0;
      0 0 0   1];

M12 = pinv(M01)*M2;
M23 = eye(4);

% Define the screw axes of the joints
S1 = [0 0 1 0 0 0]';
S2 = [0 0 1 0 -L1 0]';

% Define the twists of the joints expressed in frame {i}
A1 = Adjoint(pinv(M01)) * S1;
A2 = Adjoint(pinv(M2)) * S2;

% Define the gravity vector
g = [0 -9.81 0]; % [m/s^2]
%g = [0 0 0];

% Define the spatial inertia matrix for the two links
G1 = [zeros(3,6);
      zeros(3,3) m1*eye(3)];
G2 = [zeros(3,6);
      zeros(3,3) m2*eye(3)];
  
% Define the velocity of each link
% Ground:
V0 = [0 0 0 0 0 0]';
V0dot = [0 0 0 0 g 0]';

% Link 1:
V1 = [0 0 0 0 0 0]';
V1dot = [0 0 0 0 0 0]';

% Link 2:
V2 = [0 0 0 0 0 0]';
V2dot = [0 0 0 0 0 0]';

% Wrench applied to the end effector:
Ftip = [0 0 0 0 0 0]';

% Calculate the inverse dynamics
% thetalist = [0 0]';
% dthetalist = [0 0]';
% ddthetalist = [0 0]';
Mlist = cat(3, M01, M12, M23);
Glist = cat(3, G1, G2);
Slist = [S1 S2];


% tau = InverseDynamics(thetalist, dthetalist, ddthetalist, ...
%                       g, Ftip, Mlist, Glist, Slist);

thetalist = [linspace(0,pi/2,10);
             linspace(0,pi/2,10)];

dthetalist = diff(thetalist,1,2);
dthetalist = [dthetalist dthetalist(:,end)];
ddthetalist = zeros(2,10);

Ftipmat = zeros(size(thetalist, 2), 6);

tau = InverseDynamicsTrajectory(thetalist', dthetalist', ddthetalist', ...
                                g, Ftipmat, Mlist, Glist, Slist);
                            
                  
%% Forward Dynamics
thetalist = [0 0]';
dthetalist = [0.1745 0.1745]';
%taumat = 10*ones(100,2);
%taumat = repmat(tau', [10 1]);
%Ftipmat = zeros(size(taumat, 1), 6);
dt = 0.05;
intRes = 100;
[thetamat, dthetamat] = ForwardDynamicsTrajectory(thetalist, dthetalist, tau, g, ...
                                                  Ftipmat, Mlist, Glist, Slist, dt, ...
                                                  intRes);
                                              
for ii = 1 : 10
    ii
    robot.teach(thetamat(ii,:));
end