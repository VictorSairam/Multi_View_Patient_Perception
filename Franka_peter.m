% RBE 501 - Robot Dynamics - Fall 2021
% Worcester Polytechnic Institute
% Team 10 Project Code
%
% Instructor: L. Fichera <lfichera@wpi.edu>
% Last modified: 12/13/2021

clear, clc, close all
addpath('utils');
addpath('function library');
%% Create the manipulator
mdl_panda;

panda.tool = eye(4); % take the tool off the end effector

q = zeros(1,7);
f1 = figure(1); f1.Position = [1000 200 800 600]; hold on;
axis([-1.5 1.5 -1.5 1.5 -1 1.5])
set(gca, 'XDir','reverse'); set(gca, 'YDir','reverse');
panda.plot(q,'trail','ro','jointcolor', [.9 .9 .9],'linkcolor', [.9 .9 .9],'lightpos', [20 0 0]); %display options here
%% Create Patient Body
%cylinder
R = 0.1;
L= 1.5;
ncs = 100;
nNodes = 100;

r = R*ones(1,nNodes);
th = linspace(0,2*pi,nNodes);

%polar
y = linspace(-L/2,L/2,ncs)';

[x,z] = pol2cart(th,r);

X = repmat(x,ncs,1);
Z = repmat(z,ncs,1);
Y = repmat(y,1,nNodes);

x_lid = zeros(2,nNodes);
y_lid = repmat([-L/2;L/2],1,nNodes);
z_lid = zeros(2,nNodes);

X = [x_lid(1,:); X; x_lid(2,:)]+.35;
Y = [y_lid(1,:); Y; y_lid(2,:)];
Z = [z_lid(1,:); Z; z_lid(2,:)]-.1;

s = surf(X,Y,Z,FaceColor = 'yellow');
s.EdgeColor = 'none';

%% Calculate the forward kinematics using the Product of Exponentials
% Let us calculate the screw axis for each joint
% Put all the axes into a 6xn matrix S, where n is the number of joints
L1 = .333; L2 = .316; L3 = .0825; L4 = .384; L5 = 0.088; L6 = 0.107;

S = [0 0 1 0 0 0;
     0 1 0 -L1 0 0;
     0 0 1 0 0 0;
     0 -1 0 L1+L2 0 -L3;
     0 0 1 0 0 0;
     0 -1 0 L1+L2+L4 0 0;
     0 0 -1 0 L5 0]';

% Let us also calculate the homogeneous transformation matrix M for the
% home configuration
M = [1 0 0 L5;
     0 -1 0 0;
     0 0 -1 L1+L2+L4-L6;
     0 0 0 1];
 
%Compute froward kinematics for given S,M and q using fkine function
T = fkine(S,M,q);

%% Test foward kinematics
nTests = 20;
plotOn = false;

fprintf('---------------------Forward Kinematics Test---------------------\n');
fprintf(['Testing ' num2str(nTests) ' random configurations.\n']);
fprintf('Progress: ');
nbytes = fprintf('0%%'); 
qlim = panda.qlim; 
% Test the forward kinematics for 100 random sets of joint variables
for ii = 1 : nTests
    fprintf(repmat('\b',1,nbytes));
    nbytes = fprintf('%0.f%%', ceil(ii/nTests*100));
    
    % Generate a random configuration
    q = [qlim(1,1) + (qlim(1,2) - qlim(1,1)) * rand(), ...
         qlim(2,1) + (qlim(2,2) - qlim(2,1)) * rand(), ...
         qlim(3,1) + (qlim(3,2) - qlim(3,1)) * rand(), ...
         qlim(4,1) + (qlim(4,2) - qlim(4,1)) * rand(), ...
         qlim(5,1) + (qlim(5,2) - qlim(5,1)) * rand(), ...
         qlim(6,1) + (qlim(6,2) - qlim(6,1)) * rand(), ...
         qlim(7,1) + (qlim(7,2) - qlim(7,1)) * rand()];
    
    % Calculate the forward kinematics
    T_test = fkine(S,M,q);
    
    if plotOn
        panda.teach(q);
        title('Forward Kinematics Test');
    end
    
    assert(all(all(abs(double(panda.fkine(q)) - T_test) < 1e-10)));
end
 
fprintf('\nTest passed successfully.\n');

%% Build Path Target Poses
q = zeros(1,7);
panda.teach(q); hold On;

% Generate and display the path that the robot has to trace
t1 = 20;
t2 = 160;
a= 0.25;   %circle amplitude
d = .35;   %patient offset x
w = .4;  %length of scan
h = -.1;  %patient height
phi = pi/2;

n1 = 40;
t = linspace(t1, t2, n1)*pi/180;
x1 = d+a*(-cos(t));
y1 = w/2*ones(1,n1);
z1 = a*sin(t)+h;

R1 = zeros(3,3,n1);
for i = 1:n1
    R1(:,:,i) = axisangle2rot([0 1 0],t(i)+phi);
end

n2 = 20;
x2 = x1(length(x1))*ones(1, n2);
y2 = linspace(w/2,-w/2,n2);
z2 = z1(length(z1))*ones(1, n2);
R2 = repmat(R1(:,:,length(R1)),1,1,n2);

t = linspace(t2, t1, n1)*pi/180;
x3 = d+a*(-cos(t));
y3 = -w/2*ones(1,n1);
z3 = a*sin(t)+h;

R3 = zeros(3,3,n1);
for i = 1:n1
    R3(:,:,i) = axisangle2rot([0 1 0],t(i)+phi);
end
path = [[x1 x2 x3]; [y1 y2 y3]; [z1 z2 z3]];
pose = cat(3,R1,R2,R3);

scatter3(path(1,:), path(2,:), path(3,:), 'filled');

% Convert Cartesian coordinates into twists
targetPose = zeros(6,size(path,2)); % each column of this matrix is a target pose represented by a twist
targetHT = repmat(zeros(4,4),1,1,size(path,2));
for ii = 1 : size(path,2)
    % First calculate the homogeneous transformation matrix representing
    % each target pose
    R = pose(:,:,ii);
    T = [R path(:,ii); 
         0 0 0 1];
     
    % Then perform the matrix logarithm operation to convert transformation
    % matrices into 4x4 elements of se(3)
    t = MatrixLog6(T);
    
    % Finally, "unpack" this matrix (i.e., perform the inverse of the
    % bracket operator)
    targetHT(:,:,ii) = T;
    targetPose(:,ii) = [t(3,2) t(1,3) t(2,1) t(1:3,4)']';
end

fprintf('Built Path Target Poses\n')

%% Calculate Path Inverse Kinematics to Obtin Path Joint Variables

% Calculate the twist representing the robot's home pose
currentPose = MatrixLog6(M);
currentPose = [currentPose(3,2) currentPose(1,3) currentPose(2,1) currentPose(1:3,4)']';

% Set the current joint variables to home configuration
currentQ = [0,0,0,0,0,0];

qList = zeros(7,width(targetPose));
%Calculate inverse kinematics for each target pose
for target = 1:width(targetPose)
    
    %Set target twist pose
    V_desired = targetPose(:,target);
    
    %use petercork's ikin algorythm
    currentQ= IKinSpace(S, M,targetHT(:,:,target), currentQ(:,1), 1e-3, 1e-3);

    qList(:,target) = mod(currentQ(:,1),2*pi);
    
end

fprintf('\nTarget Joint Variables Calculated');

%% Trace Path Trajectory Using Joint Variables list

fprintf('\nAnimating Robot Path Trace.....')
%panda.plot(qList', 'trail', {'r', 'LineWidth', 5});
fprintf('\nPath Trace Completed\n')

%% Dynamics Integration Preliminaries

%masses
m(1) = 0.629769;  % [kg] Mass of Link 1
m(2)=  4.970684;  % [kg] Mass of Link 2
m(3) = 0.646926; % [kg] Mass of Link 3
m(4) = 3.228604; % [kg] Mass of Link 4
m(5) = 3.587895; % [kg] Mass of Link 5
m(6) = 1.225946; % [kg] Mass of Link 6
m(7) = 1.666555; % [kg] Mass of Link 7
m(8)= 7.35522e-01+.45; % [kg] Mass of Link 8 (adding .45 to account for endeffector)

r= 0.05; %cylinder radius in meters
%cylinder heights in meters
h(1) = .111;
h(2) = .222;
h(3) = .166;
h(4) = .15;
h(5) = .1;
h(6) = .284;
h(7) = .11;
h(8) = .03;

%calculate spatial inertia matrixes
G = zeros(6,6,8);
Z = zeros(3,3);
for i = 1:8
    I(1) = m(i)*(3*(r^2)+(h(i)^2))/12;
    I(2) = m(i)*(3*(r^2)+(h(i)^2))/12;
    I(3) = m(i)*(r^2)/2;
    Ir = diag(I);
    G(:,:,i) = [Ir Z; Z m(i)*eye(3)];
    G(:,:,i) = [Ir Z; Z m(i)*eye(3)];
end

%home configuration link frames
M01 = [1 0 0 0; 0 1 0 0; 0 0 1 h(1)/2; 0 0 0 1];
M12 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(1)+h(2))/2; 0 0 0 1];
M23 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(2)+h(3))/2; 0 0 0 1];
M34 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(3)+h(4))/2; 0 0 0 1];
M45 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(4)+h(5))/2; 0 0 0 1];
M56 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(5)+h(6))/2; 0 0 0 1];
M67 = [1 0 0 0.088; 0 -1 0 0; 0 0 -1 0; 0 0 0 1];
M78 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(7)+h(8))/2; 0 0 0 1];
M89 = [1 0 0 0; 0 1 0 0; 0 0 1 (h(8))/2; 0 0 0 1];

Mlist = cat(3, M01, M12, M23, M34, M45, M56, M67, M78, M89);
%% Dynamics
fprintf('\nCalculating Dynamics.....')
% Initialize the arrays where we will accumulate the output of the robot
% dynamics, so that we can display it later
qtt = []; % Joint Variables
qdtt =[]; %joint velocities
tau = []; %joint torques
nPts = size(qList,2);
n = size(qList,1);
g = [0 0 -9.81]; % Gravity Vector [m/s^2]

for jj = 1 : nPts - 1
    t0 = 0; tf = .5; % Starting and ending time of each trajectory
    N = 50;          % Number of intermediate setpoints
    t = linspace(t0, tf, N); % time vector
    
    q = zeros(n,N);   % joint variables
    qd = zeros(n,N);  % joint velocities
    qdd = zeros(n,N); % joint accelerations
    for ii = 1 : n
        % Calculate the coefficients of the quintic polynomial
        a = quinticpoly(t0, tf, qList(ii,jj), qList(ii,jj+1),0, 0, 0, 0);
        
        % Generate the joint profiles (position, velocity, and
        % acceleration)
        q(ii,:) = a(1) + a(2) * t + a(3) * t.^2 + a(4) * t.^3 + a(5) * t.^4 + a(6) * t.^5;
        qd(ii,:) = a(2) + 2*a(3)*t + 3*a(4)*t.^2 + 4*a(5)*t.^3 + 5*a(6)*t.^4;
        qdd(ii,:) = 2*a(3) + 6*a(4)*t + 12*a(5)*t.^2 + 20*a(6)*t.^3;
    end
    
    % Use the equations of motion to calculate the necessary torques to trace
    % the trajectory
    Ftipmat = zeros(N,6); % no end effector force
    taumat = InverseDynamicsTrajectory(q', qd', qdd', ...
        g, Ftipmat, Mlist, G, S);
    
    % Use the Forward Dynamics to simulate the robot behavior
    dt = tf/N;  % time step
    intRes = 1; % Euler integration constant
    [qt, qdt] = ForwardDynamicsTrajectory(q(:,1), qd(:,1), taumat, g, ...
        Ftipmat, Mlist, G, S, dt, ...
        intRes);
    
    qtt = [qtt; qt]; % Accumulate the results
    qdtt = [qdtt; qdt];
    tau = [tau; taumat];
end

fprintf('Done.\n');
%% 

fprintf('Simulate the robot...');
title('Inverse Dynamics Control');
panda.plot(qtt(1:10:end,:),'jointdiam', 1);
fprintf('Done.\n');
%% 

% Display the Joint Torques
color_list = [[0, 0.4470, 0.7410];	  
          	[0.8500, 0.3250, 0.0980];	         
          	[0.9290, 0.6940, 0.1250];	        
          	[0.4940, 0.1840, 0.5560];	          
          	[0.4660, 0.6740, 0.1880];	          	
          	[0.3010, 0.7450, 0.9330];
            [0, 0.5, 0]];
for kk = 1:n
    %torques
    figure(2)
    subplot(n,1,kk)
    plot((1:length(tau))*dt, tau(:,kk), 'Linewidth', 2,'color', color_list(kk,:)); 
    grid on
    xlim([0 max((1:length(tau))*dt)]);
    xlabel('Time [s]'), ylabel('Torque [Nm]');
    title_str = append('Torque Profile for Joint ',string(kk));
    title(title_str);
    %positions
    figure(3)
    subplot(n,1,kk)
    plot((1:length(qtt))*dt, qtt(:,kk), 'Linewidth', 2,'color', color_list(kk,:)); 
    grid on
    xlim([0 max((1:length(qtt))*dt)]);
    xlabel('Time [s]'), ylabel('Joint Position [radians]');
    title_str = append('Position Profile for Joint ',string(kk));
    title(title_str);
    %velocity
    figure(4)
    subplot(n,1,kk)
    plot((1:length(qdtt))*dt, qdtt(:,kk), 'Linewidth', 2,'color', color_list(kk,:)); 
    grid on
    xlim([0 max((1:length(qdtt))*dt)]);
    xlabel('Time [s]'), ylabel('Angular Velocity [radians/s]');
    title_str = append('Velocity Profile for Joint ',string(kk));
    title(title_str);
end

fprintf('Program completed successfully.\n');





