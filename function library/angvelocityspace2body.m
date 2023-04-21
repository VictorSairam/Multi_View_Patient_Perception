function omega_b = angvelocityspace2body(omega_s,R)
%angvelocityspace2body converts anvular velocity vector to body frame
% omega_b = angvelocityspace2body(omega_s,R)
    omega_b = R'*omega_s;
end