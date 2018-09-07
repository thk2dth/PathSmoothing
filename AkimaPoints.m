function [pts, der1, der2, cur] = AkimaPoints( knotVector, controlPoints, u )
% Calculate points on an Akima curve.
% input:
%   knotVector, knot vector of the Akima spline curve.
%   controlPoints, control points of the Akima spline curve.
%       number of control points equals is four times of the number of the
%       knotVector - 1. Each sub Akima curve has foour control points.
%   u, curve parameters where points are derived.
% output:
%   pts, points at curve parameters. Each column is a point.
%   der1, der2, first and second orders of derivative.
%   cur, curvature.
% HJ, 20180119
num = length(u);
dim = size(controlPoints, 1);
if nargout > 0
    pts = zeros(dim, num);
end
if nargout > 1
    der1 = zeros(dim, num);
end
if nargout > 2
    der2 = zeros(dim, num);
end
if nargout > 3
    cur = zeros(1, num);
end

for i = 1 : num
    index = findAkimaKnotVector(u(i), knotVector);
    ic = (index - 1) * 4; % Each sub curve has four control points.
    du = u(i) - knotVector(index);
    if nargout > 0
        pts(:, i) = controlPoints(:, ic+1:ic+4) * [1; du; du*du; du*du*du];
    end
    if nargout > 1
        der1(:, i) = controlPoints(:, ic+1:ic+4) * [0; 1; 2*du; 3*du*du];
    end
    if nargout > 2
        der2(:, i) = controlPoints(:, ic+1:ic+4) * [0; 0; 2; 6*du];
    end
    if nargout > 3
        cur(i) = Curvature(der1(:, i), der2(:, 2) );
    end
end
end

function index = findAkimaKnotVector(u, knots)
% Use binary search method to find the index where knots[index]<= u <
% knots[index+1].
low = 1;
high = length(knots);
mid = floor( 0.5*(high + low) );
while (high - low) > 1
    if knots(mid) > u
        high = mid;
    else
        low = mid;
    end
    mid = floor( 0.5*(low + high) );
end
index = mid;
end