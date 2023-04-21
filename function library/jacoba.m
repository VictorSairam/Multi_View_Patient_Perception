function J_a = jacoba(S,M,q)   
%JACOBE computes the analytical jacobian.
%
%   J_a = jacoba(S,M,q) 

%   Where,
%   S is a 6xn matrix containing the screw axis of the n-dof robot as its
%   colounms
%   
%   M is the 4x4 home configuration transformation
%
%   q is the nx1 vector containing the joint variables
%
%   J_a is the 3xn matrix giving the body jacobian of the n-dof robot

    J_b = jacobe(S,M,q);
    T = fkine(S,M,q);
    R = T(1:3,1:3);
    J_bv = J_b(4:6,:);
    J_a = R*J_bv;
end