function controlPoints = BsplineFittingPIA(rawData, degree, knotVector, curveParameter, tol, maxItr)
% curve fitting with PIA method
% input:
%   rawData (dim * number), data points to fit
%   degree, degree of the B-spline curve
%   knotVector, knot vector of the B-spline basis.
%   curveParameter, parameters corresponding to the data points.
%   tol, fitting tolerance.
%   maxItr, maximum number for iteration.
% output:
%   controlPoints, control points of the final bspline.
% HJ, 20160607
validateattributes(knotVector, {'numeric'}, {'real', 'vector', 'nondecreasing'});
validateattributes(curveParameter, {'numeric'}, {'real', 'vector', 'nondecreasing'});
if nargin < 5
    tol = 0.1; % default fitting tolerance.
end
if nargin < 6
    maxItr = 100; % default maximum number for iteration.
end
itrNum = 1;
basisMatrix = MatrixBspline( curveParameter, degree, knotVector );
omega = 2 / ( 1+min(eig(basisMatrix))); % omega = 2 / (1+lamda_n(B))
controlPoints = rawData;
while(1)
    pts = bspeval(degree, controlPoints, knotVector, curveParameter);
    delta = rawData - pts;
    maxErr = max( sqrt(sum(delta.^2, 1))); % maximum fitting error
    if maxErr < tol, break; end
    controlPoints = controlPoints + omega*delta;
    itrNum = itrNum + 1;
    if itrNum > maxItr, break; end
end
end