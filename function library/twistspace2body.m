function V_b = twistspace2body(V_s,T)
%twist2body converts space frame twist, V_s to V_b in the body frame .
%   V_b = twistspace2body(V_s,T)
%   
    R = T(1:3,1:3);
    P = T(1:3,4);
    Ps = skew(P);
    Ad = [R zeros(3); Ps*R R];
    V_b = pinv(Ad)*V_s;
        function S = skew(a)
            S = [    0 -a(3)  a(2);
                  a(3)     0 -a(1);
                 -a(2)  a(1)     0];
        end
end