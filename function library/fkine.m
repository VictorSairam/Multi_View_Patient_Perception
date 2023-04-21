function T = fkine(S,M,q,frame)
%FKINE Computes foward kinematics using the product of exponentials.
%   T = fkine(S,M,q,frame)
%
%   Where; 
%   T is the 4x4 homogeneus transfrom matrix from the base frame to
%   the end effector
%
%   S is a 6xn matrix containing the screw axis of the n-dof robot
%
%   M is a 4x4 homogeneus transfrom matrix from the base frame to the end
%   effector when the bot is in the home configuration corresponding to all
%   joint variables = zero
%
%   q is the nx1 vector containing the joint variables
%   takes optional variable frame
%   frame is either 'space' or 'body'
%   deafults to space frame if not specified
%
%   See also TWIST2HT, AXISANGLE2ROT, SKEW

    N = width(S);                       %Finds the quantity of joints
    Tm=1;                               %Initilizes Tm, which will be the recursive product
    
%Multiply transfrom matrixes from each joint frame together
    for i = 1:N                        
        Ti = twist2ht(S(:,i),q(i));     %Computes Homogeneus transfrom from frame i-1 to i
        Tm=Tm*Ti;                       %Multiplies new transfrom matrix Ti to the recursive product
    end
    if ~exist('frame','var') || strcmp(frame,'space')
              T=Tm*M;       %Multiplies the recursive product by the home configuration matrix
    elseif strcmp(frame,'body')
              T=M*Tm;  
    end               
end
