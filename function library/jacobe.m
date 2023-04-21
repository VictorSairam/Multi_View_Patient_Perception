function J_b = jacobe(S,M,q) 
%JACOBE computes the body jacobian.
%
%   J_b = jacobe(S,M,q) 

%   Where,
%   S is a 6xn matrix containing the screw axis of the n-dof robot as its
%   colounms
%   
%   M is the 4x4 home configuration transformation
%
%   q is the nx1 vector containing the joint variables
%
%   J_b is the 6xn matrix giving the body jacobian of the n-dof robot

     J = jacob0(S,q);
     T = fkine(S,M,q);
     J_b = zeros(6,width(S));
     for i = 1:width(J)
         J_b(:,i) = adjoint(J(:,i),pinv(T));
     end
end