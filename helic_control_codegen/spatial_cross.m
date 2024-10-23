function Vx = spatial_cross(V)

v = V(1:3);
w = V(4:6);
Vx = zeros(6,6);
skw = skew(w);
Vx(1:3,1:3) = skw;
Vx(4:6,4:6) = skw;
Vx(1:3,4:6) = skew(v);
end