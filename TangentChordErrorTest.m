function [flag, d] = TangentChordErrorTest( p0, p1, t0, t1, ce )
% Tangent-chord error test for three points.
% Tangent-chord error test is well described in Wang2014IJMTM.
% Note that this test can only be applied to planer situations.
% Input:
%   p0, p1, two points. Column-vectors.
%   t0, t1, two unit vectors. Column-vectors.
%   ce, chord error tolerance.
% Output:
%   flag, flag indicates whether the three points meet the test or not.
%       0, no; 1, yes.
%   d1, d2, chord errors for the first and second arcs.
% HJ, 20180122.
v = p1 - p0;
A = [t0, -t1]; % p0+t0*x0 = p1+t1*x1.
x = A \ v; % For planar situation, it is the same as Eq. (9) in Wang2014IJMTM.
m = p0 + t0 * x(1); % point M
d = distance2Line(p0, p1, m);
if d > ce
    flag = 0;
else
    flag = 1;
end

end

