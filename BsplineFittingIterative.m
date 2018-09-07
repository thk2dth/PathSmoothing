function [ctrl, knots] = BsplineFittingIterative( rawData, p, ce, N )
% Bspline fitting with chord error constraints.
% First, interpolate the data points. Then, for the curve segment
% [u_i, u_i{i+1}] that violates the chord error constraint, insert new knot
% vector 0.5*(u_i+u_i+1). The points corresponding to the inserted knot
% is middle point of p_i and p_i+1.
% Input:
%   rawData, data to be fitted. Each column is a point.
%   p, degree of the B-spline curve.
%   ce, chord error tolerance.
%   N, number of points on the B-spline curve used to check the chord
%       error.
% Output:
%   ctrl, control points.
%   knots, knot vector.

% HJ, 20180122.

[~, n] = size(rawData);
% [knots, cp] = CalculateKnotVector(rawData, p); % knot vector and curve parameters
% knotsNew = knots;
dataNew = rawData;
% cpNew = cp;
numItr = 0;
numMaxItr = 10; % maximum number of iteration.
while(numItr <= numMaxItr)
    [knots, cp] = CalculateKnotVector(rawData, p); % knot vector and curve parameters
    bm = MatrixBspline( cp, p, knots ); % basis matrix,
    % bm is size number of parameters * number of data points;
    % rawData is size dimension * number of data points.
    ctrl = rawData / bm'; % control points.
    ni = 0; % number of insert operations.
    for i = 1 : n-1
        u = linspace(cp(i), cp(i+1), N);
        tp = bspeval(p, ctrl, knots, u);
        for j = 1:N
            d = Distance2Line(rawData(:, i), rawData(:, i+1), tp(:, j) );
            if d > ce
                %                 uAdd = 0.5 * (cp(i) + cp(i+1) ); % insert the middle curve parameter as the knot.
                dataAdd = (rawData(:, i) + rawData(:, i+1) ) * 0.5;
                % index of the inserted knot.
                %                 ik = findspan(n+ni-1, p, uAdd, knotsNew) + 1; % refer to function findspan()
                %                 knotsNew = [knotsNew(1:ik), uAdd, knotsNew(ik+1 : n+ni+p+1)];
                dataNew = [dataNew(:, 1:i+ni), dataAdd, dataNew(:, i+ni+1:n+ni)];
                %                 cpNew = [cpNew(1:i+ni), uAdd, cpNew(i+ni+1:n+ni)];
                ni = ni + 1;
                break;
            end
        end
    end
    if ni == 0
        return;
    else
        %         cp = cpNew;
        %         knots = knotsNew;
        rawData = dataNew;
        numItr = numItr +1;
        n = n + ni; % number of data needs update too.
    end
end
end


