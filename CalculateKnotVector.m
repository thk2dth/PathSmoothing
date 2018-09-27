function [ knotVector, curveParameter ] = CalculateKnotVector( rawData, p)
% calculate the curve parameter corresponding to the data points by chord length firstly;
% then the knot vector is obtained.
%  see Chapter 9.2 in "The NURBS Book" . 
% input:
%   rawData (dim * number), data points to fit
%   p, degree of the B-spline curve
% output:
%   knotVector, knot vector of the B-spline
%   curveParameter, parameters corresponding to the data points.
% HJ, 20160606
validateattributes(rawData, {'numeric'}, {'real'});
validateattributes(p, {'numeric'}, {'positive','integer','scalar'});
number = size(rawData, 2);
rawDataDiff = diff(rawData, 1, 2);
chordLength = sqrt( sum( rawDataDiff.^2, 1));
curveParameter = zeros(1, number);
curveParameter(number) = 1.0; % suppose the knot vector interval is [0, 1].
knotVector = zeros(1, number+p+1);
knotVector(number+1:end) = 1.0;
lengthSum = sum(chordLength);
for i=2:number-1
    curveParameter(i) = curveParameter(i-1) + chordLength(i-1)/lengthSum;
end
for i =1:number-p-1
    knotVector(i+p+1) = sum(curveParameter(i+1:i+p)) / p;
end
end

