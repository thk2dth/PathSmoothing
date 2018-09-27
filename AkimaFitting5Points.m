function [ controlPoints, knotVector ] = AkimaFitting5Points( rawData )
% Fit five points with an Akima spline.
% Follow the description in Wang et. al., A real-time look-ahead interpolation algorithm
% based on Akima curve fitting, IJMTM, 2014.
% input:
%   rawData, (dim * number), data points to fit. Each column is a point.
% output:
%   knotVector, knot vector of the Akima spline curve.
%   controlPoints, control points of the Akima spline curve.
%       number of control points is four times of the number of the
%       points - 1. Each sub Akima curve has four control points.
% HJ, 20180119
[dim, num] = size(rawData); % dimension and number of data points.
knotVector = CalculateKnots(rawData);
controlPoints = zeros(dim, (num-1) * 4 );
U = zeros(dim, num-1+4); % total number U vector is num-1+4.
for i = 1 : num-1 % Eq. (3).
    U(:, i+2) = ( rawData(:, i+1) - rawData(:, i) ) / ...
        ( knotVector(i+1) - knotVector(i) );
end
U(:, 2) = 2 * U(:, 3) - U(:, 4); % U_{-1}
U(:, 1) = 2 * U(:, 2) - U(:, 3); % U_{-2}
U(:, num+2) = 2 * U(:, num+1) - U(:, num); % U_n
U(:, num+3) = 2 * U(:, num+2) - U(:, num+1); % U_{n+1}
for i = 1:num - 1 % There are num-1 Akima curves.
    ic = 4 * (i-1);
    ui = U(:, i+2); % current slope
    controlPoints(:, ic + 1) = rawData(:, i); % A = P_i
    s1 = ApproximateDer( U(:, i:i+3) );
    controlPoints(:, ic + 2) = s1; % B = s(u_i)
    s2 = ApproximateDer( U(:, i+1:i+4) );
    du = knotVector(i+1) - knotVector(i); % u_{i+1} - u_{i}
    controlPoints(:, ic + 3) = (ui*3 - s1*2 - s2 ) / du; % C
    controlPoints(:, ic + 4) = (-ui*2 + s1 + s2) / (du * du); % D
end

end

% Calculate the knots according to centripetal parameterization method.
function knots = CalculateKnots( rawData )
num = size(rawData, 2); % number of data points.
knots = zeros(1, num); % number of knots equals to the number of data points.
p = 4; % p-norm is used to calculate the length.
for i = 2 : num
    v = rawData(:, i) - rawData(:, i-1);
    knots(i) = norm(v, p) + knots(i-1); % accumulative length.
end
knots = knots / knots(num); % normalize the knots.
end

% Approximate the derivative.
function si = ApproximateDer( U )
l1 = norm( U(:, 2) - U(:, 1) );
l2 = norm( U(:, 4) - U(:, 3) );
si = ( U(:, 2)*l2 + U(:, 3)*l1 ) / (l1 + l2); 
end