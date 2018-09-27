function [ ctrl, knots, degree ] = BsplineTransition( rawData, ce, r )
% Bspline transition for linear tool path.
% After transition, the remaining linear tool path and the B-spline
% curves will be merged together as a united B-spline.
% Input:
%   rawData, data points.
%   ce, chord error.
%   r, d1 / d2. r = 0.25, by default.
% Output:
%   ctrl, control points of the united B-spline.
%   knots, knot vector of the united B-spline, belonging to [0, 1].
%   degree, degree of the united B-spline. 3, by default.
% HJ, 20180125.

% The inserted B-spline has five control points, and the B-spline
% converted by the remaining linear tool path has four control points.
% Therefore, for the united B-spline, there will be 4*4+5*3-6=25 control
% points, and 8*4+9*3-6*5=29 knots.
degree = 3;
[dim, num] = size(rawData);
ctrl = zeros(dim, 25);
knots = zeros(1, 29);
d = zeros(1, 3); % three d2 values, for inserting the intermediate points.
% First, place the knots and control points of the inserted B-spline in the
% final knots and control points.
for i = 2 : num-1
    [c, k, d2] = Transition3Points(rawData(:, i-1), rawData(:, i), ...
        rawData(:, i+1), ce, r);
    % There are total 7 segments, 3 B-splines and 4 remaining linear
    % segments.
    % the multiplicity of the parameter at the junction is 3, not 4.
    knots(7*i-9 : 7*i-6) = (k(2:5) + (2*i-3) ) / 7.0;
    ctrl(:, 7*i-10 : 7*i-6) = c(:, 1:5); % the first four points.
    d(i-1) = d2;
end
% Then, place the knots and control points of the converted B-spline in the
% final knots and control points.
% There are four remaining segments.
for i = 1 : num-1
    if 1 == i
        ps = rawData(:, 1); % the start point of the first remaining segment.
        d1p = 0;
    else
        ps = ctrl(:, 7*i-6);
        d1p = d(i-1);
    end
    if num-1 == i
        pe = rawData(:, num); % the end point of the last remaining segment.
        d2p = 0;
    else
        pe = ctrl(:, 7*i-3);
        d2p = d(i);
    end
    vp = pe - ps;
    lp = sqrt(vp' * vp);
    % Since there is no length margin used during transition, we do not try
    % to make the concatenated curve C^1 continuous.
    % The two intermediate points are just even distributed in the remaining
    % line segment.
    v1 = ps + vp * (lp / 3);
    v2 = pe - vp * (lp / 3);
    
    if 1 == i
        ctrl(:, 1 : 3) = [ps, v1, v2];
    elseif num-1 == i
        ctrl(:, 23 : 25) = [v1, v2, pe];
    else
        ctrl(:, 7*i-5 : 7*i-4) = [v1, v2];
    end
    if 1 == i
        knots(1:4) = [0, 0, 0, 0];
    elseif num-1 == i
        knots(23:29) = [6, 6, 6, 7, 7, 7, 7] / 7.0;
    else
        knots(7*i-5:7*i-3) =[i-1, i-1, i-1] * (2.0/ 7.0);
    end
end

end

function [c, k, d2] = Transition3Points(p1, p2, p3, ce, r)
v1 = p2 - p1;
v2 = p3 - p2;
l1 = sqrt(v1' * v1);
l2 = sqrt(v2' * v2);
m1 = v1 / l1;
m2 = v2 / l2;
beta = 0.5 * acos( m1' * m2 );
% Eq. (9)
n1 = 2 * ce / sin(beta);
n2 = min(l1, l2) / (2*(1+r) ); % lm = min(l1, l2)/2.
d2 = min(n1, n2);
% inserted B-spline has five control points.
c = zeros(3, 5);
k = [0, 0, 0, 0, 0.5, 1, 1, 1, 1];
c(:, 1) = p2 - m1 * (d2*(1+r) );
c(:, 2) = p2 - m1 * d2;
c(:, 3) = p2;
c(:, 4) = p2 + m2 * d2;
c(:, 5) = p2 + m2 * (d2*(1+r) );
end