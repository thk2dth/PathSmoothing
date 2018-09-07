function [flag, d1, d2] = BiChordErrorTest( p0, p1, p2, ce )
% Bi-chord error test for three points.
% Bi-chord error test is well described in Wang2014IJMTM, and
% Yau2007IJMTM.
% Input:
%   p0, p1, p2, three points. Column-vectors.
%   ce, chord error tolerance.
% Output:
%   flag, flag indicates whether the three points meet the test or not.
%       0, no; 1, yes.
%   d1, d2, chord errors for the first and second arcs.
% HJ, 20180122.
v1 = p1 - p0;
v2 = p2 - p1;
l1 = sqrt(v1' * v1);
l2 = sqrt(v2' * v2);
% acos return values belong to [0, pi].
tc = v1' * v2 / (l1*l2); % temporary variable.
alpha = acos(tc); % alpha = pi- theta.
phi = atan(l1 * sin(alpha) / (l2 + l1*tc) );
if phi < 0
    phi = phi + pi; % if phi<0, the arc determined by p0 and p1 is a major arc.
end
% r is always greater than 0.
r = 0.5 * l1 / sin(phi);
d1 = r * (1 - cos(phi) );
d2 = r * (1 - cos(alpha - phi) );
if d1 > ce || d2 > ce
    flag = 0;
else
    flag = 1;
end

end

