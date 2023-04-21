function J = jacob0(S,q)
%JACOB0 Computes Space Jacobian.
%   J = jacob0(S,q)
%   
%   Where,
%   S is a 6xn matrix containing the screw axis of the n-dof robot as its
%   colounms
%   
%   q is the nx1 vector containing the joint variables
%
%   J is the 6xn matrix giving the space jacobian of the n-dof robot
%
%   See also ADJOINT, TWIST2HT

    N = width(S);                            %Finds the quantity of joints
    Ta = zeros(4,4,N);                       %Initilizes Ta 4x4xN 3D array
    J = zeros(6,N);                          %Initilizes J 6xN space jacobian
    
%For each joint calculate its respective jacobian column
    for i=1:N
    %Find Transformation Matrix corresponding to twist i    
        Ta(:,:,i) = twist2ht(S(:,i),q(i));   
        T = eye(4);   
    %Recursivly multiply each transfrom matrix from j = 1 to i-1    
        for j = 1:i-1
            T = T*Ta(:,:,j);                 
        end
    %adjoint function transfroms twist axis i into the ith column of the space jacobian        
        J(:,i) = adjoint(S(:,i),T);          
    end
end