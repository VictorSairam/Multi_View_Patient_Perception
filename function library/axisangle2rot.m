function R = axisangle2rot(omega,theta)
%AXISANGLE2ROT Computes the rotation matrix R for given axis and angle
%   R = axisangle2rot(omega,theta)
%
%   See also FKINE, TWIST2HT, SKEW

    w = skew(omega);
%Apply Rodrigues' rotation formula
    R = eye(3,3)+sin(theta)*w + (1-cos(theta))*(w^2);
end