function S = skew(a)
%SKEW computes skew matrix, S, for the 3x1 vextor a
%   S = skew(a)

    S = [    0 -a(3)  a(2);
          a(3)     0 -a(1);
         -a(2)  a(1)     0];
end