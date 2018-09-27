function [nrb, curveParameter] = BsplineFitting( rawData, degree, method, tol )
% b-spline curve fitting with different method
% for interpolation method see "The NURBS Book"
% for PIA method see "Weighted progressive iteration approximation and convergence analysis"
% input:
%   rawData (dim * number), data points to fit
%   degree, degree of the B-spline curve
%   method, fitting
%   tol, tolerance
% output:
%   nrb, b-spline (nurbs) curve after fitting. see function nrbmak
%   curveParameter, parameters corresponding to the data points.
% HJ, 20160607
validateattributes(rawData, {'numeric'}, {'real'});
validateattributes(degree, {'numeric'}, {'positive','integer','scalar'});
if nargin > 2
    validateattributes(method, {'numeric'}, {'positive','integer','scalar'});
else
    method = 1; % use the interpolation method by default
end
[ knotVector, curveParameter ] = CalculateKnotVector( rawData, degree);
switch method
    case {1} % 'interpolation'
        basisMatrix = MatrixBspline( curveParameter, degree, knotVector );
        controlPoints = basisMatrix \ rawData'; % with size of  (number * dm)
        nrb = nrbmak(controlPoints', knotVector);
    case {2} % 'PIA'
        nrb = BsplineFittingPIA(rawData, degree, knotVector, curveParameter, tol);
    case {3} % fast and accurate fitting method.
        nrb = BsplineFittingFast(rawData, degree, knotVector, curveParameter, tol*0.2);
        
end

end