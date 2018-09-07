function [pts, der1, der2, cur] = BezierFastPoints( knotInterval, controlPoints, N )
% Calculate N points on the piecewise Bezier curve.
% input:
%   knotVector, knot vector of the Akima spline curve.
%   controlPoints, control points of the Akima spline curve.
%       number of control points equals is four times of the number of the
%       knotVector - 1. Each sub Akima curve has foour control points.
%   N, number of points to be evaluated.
% output:
%   pts, points at curve parameters. Each column is a point.
%   der1, der2, first and second orders of derivative.
%   cur, curvature.
% HJ, 20180122
p = 3; % degree of the Bezier curve.
[dim, ns] = size(controlPoints);
ns = ns/(p+1) + 2; % number of curve segments.
nt = floor(N/ns);
nps = [ones(1, ns-1)*nt, N - nt*(ns-1)]; % number of points of each segments.
npa = 1; % accumulative number of points.
if nargout > 0
    pts = zeros(dim, N);
%     B = @(t)bernsteinMatrix(p, t)'; % bernstein matrix.
    B = @(t) [-(t - 1)^3; 3*t*(t - 1)^2; -3*t^2*(t - 1); t^3];
end
if nargout > 1
    der1 = zeros(dim, N);
%     Bd1(t) = diff(B(t), 1); % first order derivative
    Bd1 = @(t) [-3*(t-1)^2; 9*t^2-12*t+3; -3*t*(3*t-2); 3*t^2];
end
if nargout > 2
    der2 = zeros(dim, N);
%     Bd2(t) = diff(B(t), 2); % second order derivative
Bd2 = @(t) [6*(1-t); 6*(3*t-2); 6*(1-3*t); 6*t];
end
if nargout > 3
    cur = zeros(1, N);
end

for i = 1 : ns
    if i == 1
        k = i + 1; % the first segment shares control points with the second segment.
    elseif i == ns
        k = i - 1; % the last segment shares control points with the previous segment.
    else
        k = i;
    end
    ctrl = controlPoints(:, (p+1)*(k-2)+1 : (p+1)*(k-1) );
    u = linspace(knotInterval(1, i), knotInterval(2, i), nps(i) );
    for j = 1 : nps(i)
        if nargout > 0
            pts(:, npa+j-1 ) = ctrl * B( u(j) );
        end
        if nargout > 1
            der1(:, npa+j-1) = ctrl * Bd1(u(j) );
        end
        if nargout > 2
            der2(:, npa+j-1) = ctrl * Bd2(u(j) );
        end
        if nargout > 3
            % absolute curvature.
            cur(npa+j-1) = Curvature(der1(:, npa+j-1), ...
                der2(:, npa+j-1) );
        end
    end
    npa = npa + nps(i);
end
end

