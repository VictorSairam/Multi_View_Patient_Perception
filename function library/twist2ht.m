function T = twist2ht(S,theta)
    omega = S(1:3);
    v  = S(4:6);
    w = skew(omega);
    R = axisangle2rot(omega,theta);
    P = ((theta*eye(3))+((1-cos(theta))*w )+((theta-sin(theta))*(w^2)))*v;
    T = [R P;  0 0 0 1]; 
end