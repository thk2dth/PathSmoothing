function [flag, der1, der2] = AlterBsplineType(pts, p, ctrl, knots, curveParameter, ce )
% Determine which one of the six middle control points should be altered.
% After knot inserting, each curve parameter in the knot vector has a multiplicity of the
% degree, which means the curve point at the parameter is one of the control
% points.
% input:
%   pts, original data points.
%   p, degree of the bspline.
%   ctrl, control points of the bspline.
%   knots, knot vector of the bspline.
%   curveParameter, the curve parameters where the curve is splitted.
%   ce, chord error.
% output:
%   flag, determine the points to be altered.
%           1 (000001), the first point; 2 (000010), the second one;
%           4 (000100), the third one; 8 (001000), the fourth one;
%           16 (010000), the fifth one; 32 (100000), the six one;
%           0 (000000), none.
%           e.g., flag = 6 implies the second and third points are out of the constraints. 
%   der1 (optional), first derivative.
%   der2 (optional), second derivative.
% HJ, 201608
% The number of intermediate points increases from 4 to 6.
%  HJ, 20180125
validateattributes(pts, {'numeric'}, {'real'});
validateattributes(curveParameter, {'numeric'}, {'real', 'vector', 'nondecreasing'});
validateattributes(ce, {'numeric'}, {'real', 'scalar'});
lenPara = length(curveParameter);
dim = size(ctrl, 1);
ns = p + 4; % number of control points per sub-curve.
flag = zeros(1, lenPara-1);
if nargout > 1
    der1 = zeros(dim, lenPara);
    if nargout > 2
        der2 = der1;
    end
end
for i = 1:lenPara-1
    for j = 2:7 % six middle points
        dis = Distance2Line(pts(:, i), pts(:, i+1),...
            ctrl(:, ns*(i-1) + j));
        if dis > ce
            flag(i) = flag(i) + 2^(j-2);
        end
    end
    if nargout > 1
        u1 = knots(ns*(i-1) + 4);
        u2 = knots(ns*(i-1) + 5);
        u3 = knots(ns*(i-1) + 6);
        der1(:, i) = p * (ctrl(:, ns*(i-1) + 2) - ctrl(:, ns*(i-1) + 1) )...
            / (u2 -  u1);
        if nargout > 2
            der2(:, i) = p*(p-1) * ( (ctrl(:, ns*(i-1) + 3) - ctrl(:, ns*(i-1) + 2)) ...
            / ((u2 - u1) * (u3 - u1)) ) - (p-1)*der1(:, i) / (u2 - u1);
        end
    end
end
% The derivatives at the last curve parameter cannot be calculated by in
% the above forward manner, and has to be derived in the backward manner.
if nargout > 1
    i = lenPara-1;
    u1 = knots(ns*(i-1) + 9);
    u2 = knots(ns*(i-1) + 8);
    u3 = knots(ns*(i-1) + 7);
    der1(:, lenPara) = p * (ctrl(:, ns*(i-1) + 8) - ctrl(:, ns*(i-1) + 7))...
        / (u1 - u2);
    if nargout > 2
        der2(:, lenPara) = -p*(p-1) * ( (ctrl(:, ns*(i-1) + 7) - ctrl(:, ns*(i-1) + 6)) ...
            / ((u1 - u2) * ( u1 - u3)) )...
            + (p-1)*der1(:, lenPara) / (u1 - u2);
    end
end
end
