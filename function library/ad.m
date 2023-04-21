function adV = ad(V)
%AD Computes the Lie Bracket of twist V;
%
%   Where v is a 6x1 twist vector
%   and adV is a 6x6 matrix

rot = V(1:3);
trans = V(4:6);

rot_skew = skew(rot);
trans_skew = skew(trans);
z = zeros(3,3);
adV = [rot_skew z; trans_skew rot_skew];

end