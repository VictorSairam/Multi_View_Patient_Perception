% Home position of the links
M01 = [eye(3) [0.05 0 0]'; 0 0 0 1];
M12 = [eye(3) [0.1 0 0]'; 0 0 0 1];
M23 = [eye(3) [0.05 0 0]'; 0 0 0 1];

M02 = M01 * M12;
M03 = M02 * M23;

% Screw Axes
S = [0 0 1 0 0 0;
     0 0 1 0 -0.1 0]';

A = zeros(6,2);
A(:,1) = Adjoint(pinv(M01)) * S(:,1);
A(:,2) = Adjoint(pinv(M02)) * S(:,2);


