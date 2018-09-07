% formulas derivation using the symbolical toolbox
% the analytical expression of the distance from a point to a line.
% HJ, 20160608.
syms x1 y1 z1 x2 y2 z2 x y z real;
pt1 = [x1; y1; z1];
pt2 = [x2; y2; z2];
pt = [x; y; z];
v1 = pt - pt1;
v2 = pt - pt2;
v3 = pt2 - pt1;
h = simplify(norm(cross(v1, v2)) / norm(v3));