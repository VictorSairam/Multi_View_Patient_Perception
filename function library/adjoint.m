function output = adjoint(V,T)
%ADJOINT Transfroms twists vectors into a new frame.
%   Vtrans = adjoint(V,T)
%   
%   Where,
%   V is the 6x1 twist in frame a
%
%   Vtrans is the 6x1 twist in frame b
%
%   T is the transfrom matrix from frame a to b
% OR
% Gives adjoint transformation %matrix for given homogeneus transform T:
%   AdT = adjoint(T)

%   See also JACOB0, SKEW
    switch nargin
        case 1
            %Get rotation and translation parts of T
            R = V(1:3,1:3);
            P = V(1:3,4);
            
            %Compute 6x6 adjoint transfromation matrix
            Ps = skew(P);
            Ad = [R zeros(3); Ps*R R];
            output = Ad;
        case 2
                %Get rotation and translation parts of T
            R = T(1:3,1:3);
            P = T(1:3,4);
            
            %Compute 6x6 adjoint transfromation matrix
            Ps = skew(P);
            Ad = [R zeros(3); Ps*R R];
            %Apply frame transformation
            Vtrans = Ad*V;
            output = Vtrans;
    end
end
