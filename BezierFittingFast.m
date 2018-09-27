function [ controlPoints, knotInterval ] = BezierFittingFast( rawData )
% Fast cubic Bezier fitting method.
% Follow the description in Yau et. al., Fast Bezier interpolator with
% real-time lookahead function for high-accuracy machining, IJMTM, 2007.
% input:
%   rawData, (dim * number), data points to fit. Each column is a point.
% output:
%   knotInterval, knot interval of each Bezier curve.
%       The number of Bezier curve segments is the number of rawData - 1.
%   controlPoints, control points of the Bezier curves.
%       The number of Bezier curves will be number of rawData - 3.
%       Each cubic Bezier curve has four control points. Thus, the 
%       total number of control points is 4 * (number of rawData - 3).
% HJ, 20180120
[dim, num] = size(rawData);
knotInterval = zeros(2, num - 1); % knot intervals of the Bezier curve segments.
p = 3; % degree of the Bezier curve.
controlPoints = zeros(dim, (p+1) * (num-3) );
vec = zeros(dim, num - 1);
len = zeros(1, num - 1);
for i = 1 : num - 1
    vec(:, i) = rawData(:, i+1) - vec(:, i);
    len(i) = norm(vec(:, i) );
end

for i = 2 : num - 2
    ic = (p+1) * (i-2);
    controlPoints(:, ic + 1) = rawData(:, i-1); % Q0
    controlPoints(:, ic + 4) = rawData(:, i+2); % Q3
    s = len(i-1) + len(i) + len(i+1);
    t1 = len(i-1) / s;
    t2 = (len(i-1) + len(i) ) / s;
    tp = bernsteinMatrix(p, [t1, t2]);
    e1 = rawData(:, i) - rawData(:, i-1) * tp(1, 1) - rawData(:, i+2) * tp(1, 4);
    e2 = rawData(:, i+1) - rawData(:, i-1) * tp(2, 1) - rawData(:, i+2) * tp(2, 4);
    tq = tp(1, 2)*tp(2, 3) -  tp(2, 2)*tp(1, 3); % temporary variable.
    controlPoints(:, ic + 2) = (e1 * tp(2, 3) - e2 * tp(1, 3) ) / tq; % Q1
    controlPoints(:, ic + 3) = (e2 * tp(1, 2) - e1 * tp(2, 2) ) / tq; % Q2
    knotInterval(:, i) = [t1; t2];
    if 2 == i
        knotInterval(:, 1) = [0; t1]; % the first Bezier curve segment.
    end
    if num-2 == i
        knotInterval(:, num-1) = [t2; 1]; % the last Bezier curve segment.
    end
end
end