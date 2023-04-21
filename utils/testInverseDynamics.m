clear; clc;

n = 2; % number of links
L1 = 0.3; % [m]
L2 = 0.4; % [m]
m1 = 1;   % [kg]
m2 = 1;   % [kg]

g = [0; -9.8; 0];

% Define the home configurations of each link, expressed in the space frame                
M1 = [eye(3) [L1 0 0]'; 0 0 0 1];
M2 = [eye(3) [L1+L2 0 0]'; 0 0 0 1];
M3 = M2;

% Define the home configurations of each link, expressed w.r.t. the
% previous link frame
M01 = M1;
M12 = [eye(3) [L2 0 0]'; 0 0 0 1];
M23 = eye(4);

Mlist = cat(3, M01, M12, M23);

% Define the screw axes of each joint, expressed in the space frame
S = zeros(6,n);
S(:,1) = [0 0 1 0 0 0]';
S(:,2) = [0 0 1 -cross([0 0 1], [L1 0 0])];

% Calculate the screw axes of each joint, expressed in the local link frame
A = zeros(6,n);
A(:,1) = Adjoint(pinv(M1)) * S(:,1);
A(:,2) = Adjoint(pinv(M2)) * S(:,2);

% % Initialize the twists and accelerations of each link
V = zeros(6,n);
Vd = zeros(6,n);
 
% Initialize the joint positions and velocities
q = zeros(2,1);
qd = ones(2,1);
qdd = ones(2,1);

%Ftip = [1; 1; 1; 1; 1; 1];
Ftip = ones(6,1);

%G1 = diag([0.010267, 0.010267, 0.00666, 3.7, 3.7, 3.7]);
%G2 = diag([0.22689, 0.22689, 0.0151074, 8.393, 8.393, 8.393]);
%G3 = diag([0.0494433, 0.0494433, 0.004095, 2.275, 2.275, 2.275]);
%Glist = cat(3, G1, G2, G3);
G1 = zeros(6,6); G1(4:6,4:6) = m1 * eye(3);
G2 = zeros(6,6); G2(4:6,4:6) = m2 * eye(3);
Glist = cat(3, G1, G2);

[~,V,Vd] = InverseDynamics(q, qd, qdd, g, Ftip, Mlist, Glist, S)



% thetalist = [0.1; 0.1; 0.1];
% dthetalist = [0.1; 0.2; 0.3];
% ddthetalist = [2; 1.5; 1];
% g = [0; 0; -9.8];
% 
% 
% M01 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.089159]; [0, 0, 0, 1]];
% M12 = [[0, 0, 1, 0.28]; [0, 1, 0, 0.13585]; [-1, 0 ,0, 0]; [0, 0, 0, 1]];
% M23 = [[1, 0, 0, 0]; [0, 1, 0, -0.1197]; [0, 0, 1, 0.395]; [0, 0, 0, 1]];
% M34 = [[1, 0, 0, 0]; [0, 1, 0, 0]; [0, 0, 1, 0.14225]; [0, 0, 0, 1]];
% 
% 
% 
% Glist = cat(3, G1, G2, G3);
% Mlist = cat(3, M01, M12, M23, M34); 
% 
% Slist = [[1; 0; 1;      0; 1;     0], ...
%        [0; 1; 0; -0.089; 0;     0], ...
%        [0; 1; 0; -0.089; 0; 0.425]];
%    
