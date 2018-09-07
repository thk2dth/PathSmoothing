function [pts, der1, der2, cur] = BsplinePoints( d, c, k, u )
% Calculate the points on a B-spline curve.
% Input:
%   d, degree of the B-spline.
%   c, control points.
%   k, knot vector.
%   u, curve parameters where the points are evaluated.
% Output:
%   pts, points at curve parameters. Each column is a point.
%   der1, der2, first and second orders of derivative.
%   cur, curvature.
% HJ, 20180125.

num = length(u); % number of curve parameters.
dim = size(c, 1);
if nargout > 0
    pts = zeros(dim, num);
end
if nargout > 1
    der1 = zeros(dim, num);
    [dc, dk] = bspderiv(d, c, k);
end
if nargout > 2
    der2 = zeros(dim, num);
    [ddc, ddk] = bspderiv(d-1, dc, dk);
end
if nargout > 3
    cur = zeros(1, num);
end

for i = 1 : num
    if nargout > 0
        pts(:, i) = bspeval(d, c, k, u(i) );
    end
    if nargout > 1
        der1(:, i) = bspeval(d-1, dc, dk, u(i) );
    end
    if nargout > 2
        der2(:, i) = bspeval(d-2, ddc, ddk, u(i) );
    end
    if nargout > 3
        cur(i) = Curvature(der1(:, i), der2(:, i) );
    end
end

end

