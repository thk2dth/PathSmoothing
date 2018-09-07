function [ctrlNew, knotsNew] = BsplineFittingFast( rawData, degree, ce, er )
% curve fitting with PIA method
% input:
%   rawData (dim * number), data points to fit
%   degree, degree of the B-spline curve
%   knotVector, knot vector of the B-spline basis.
%   curveParameter, parameters corresponding to the data points.
%   ce, chord error for curve fitting.
%   er, error ratio, error in rough fitting / tolerance.
% ouput:
%   ctrlNew, control points of the b-spline curve after fitting
%   knotsNew, knot vector of the b-spline curve after fitting.
% HJ, 20180115

validateattributes(rawData, {'numeric'}, {'real'});
validateattributes(degree, {'numeric'}, {'positive','integer','scalar'});
% Calculate the knot vector and curve parameter.
[ knots, curveParameter ] = CalculateKnotVector( rawData, degree);
% Rough fitting
ctrl = BsplineFittingPIA( rawData, degree, knots, curveParameter, ce*er );
% Bspline splitting.
[ctrlNew, knotsNew] = BsplineSplitting(degree, ctrl, knots, curveParameter);
% Determine the areas to which the local refinement will be applied.
[flag, der1, der2] = AlterBsplineType(rawData, degree, ctrlNew, knotsNew, curveParameter, ce );
ns = degree + 4; % number of control points per sub-curve.
% Local refinement.
for i = 1 : length(flag)
    if flag(i) > 0
        % knots and control points after modification.
%         fprintf('Now locally refine the %dth sub-curve...\n', i);
        % nrb.coefs is four-dimensional.
        [w, ctrl] = BsplineLocalModification( curveParameter(i), curveParameter(i+1), ...
            ctrlNew(:, (i-1)*ns+1 ), ctrlNew(1:3, i*ns+1 ), ...
            der1(:, i), der2(:, i), der1(:, i+1), der2(:, i+1), ce*(1-er) );
        ctrlNew(:, (i-1)*ns+2 :ns*i) = ctrl; % replace the altered six control points.
        knotsNew(ns*(i-1)+5 : ns*(i-1)+8) = w; % repalce the four knots.
    end  
end

