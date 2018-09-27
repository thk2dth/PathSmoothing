% Compare the consuming time of the different fitting methods.
% The comparison consists of the following steps.
% (1) Generate 5*100 points in 3-D space. The distance and angle between the adjacent
% points are both randomly generated and belong to a specified range.
% (2) The Akima fitting, the fast Bezier fitting, the iterative B-spline
% fitting, and the proposed two-step methods are applied to the generated
% 500 points. Hence there is 100 loops, and each loop processes 5 points.
% Since only the consuming time is compared, the interpolated points and
% the actual chord error are not calculated.
% HJ, 20180123.

clear all;
numItr = 100;
numPts = 5;
ce = 0.04; % chord error
er = 0.2; % error ratio. The smaller, the more stringent for rough fitting.
p = 3; % degree of B-spline.
lenRange = [0.5, 2];
alphaRange = [2.5, 30] * (pi/180); % the angle range between adjacent lines.
betaRange = [-pi, pi]; % no limit for beta.
numIterativeBspline = 50;
rawData = zeros(3, numPts, numItr); % 3-d matrix.

% Generate the 5*100 data points.
for i = 1 : numItr
    % random length of each line.
    lineLength = rand(numPts-1, 1) * (lenRange(2) - lenRange(1) ) + lenRange(1);
    alpha = rand(numPts-2, 1) * (alphaRange(2) - alphaRange(1) ) + alphaRange(1);
    % beta can be any value.
    beta = rand(numPts-2 , 1) * (betaRange(2) - betaRange(1) ) + betaRange(1);
    % the first point is at the origin.
    % the second point is random.
    t = rand(3, 1);
    rawData(:, 2, i) = t * ( lineLength(1) / sqrt(t'*t) );
    for k = 3 : numPts
        v = rawData(:, k-1, i) - rawData(:, k-2, i);
        v = v / sqrt(v'*v);
        basisPre = null(v' ); % norm vectors of v.
        rawData(:, k, i) = rawData(:, k-1, i) + [basisPre, v] * [sin(alpha(k-2) )...
            *cos(beta(k-2) ); sin(alpha(k-2) )*sin(beta(k-2) );...
            cos(alpha(k-2) )] * lineLength(k-1);
    end
end

% Apply the fitting methods to the generated data points and compare the
% consuming time.

conTime = zeros(1, 5);
% Akima fitting.
tic
for i = 1 : numItr
    % Bi-chord error test.
    for k = 2 : numPts-1
        BiChordErrorTest(rawData(:, k-1, i), rawData(:, k, i),...
            rawData(:, k+1, i), ce);
    end
    AkimaFitting5Points( rawData(:, :, i) );
end
conTime(1) = toc;

% Fast Bezier fitting.
tic
for i = 1 : numItr
    % Bi-chord error test.
    for k = 2 : numPts-1
        BiChordErrorTest(rawData(:, k-1, i), rawData(:, k, i),...
            rawData(:, k+1, i), ce);
    end
    BezierFittingFast( rawData(:, :, i) );
end
conTime(2) = toc;

% Iterative B-spline.
tic
for i = 1 : numItr
    BsplineFittingIterative( rawData(:, :, i),  p, ce,  numIterativeBspline);
end
conTime(3) = toc;

% Transition.
tic
for i = 1 : numItr
    BsplineTransition( rawData(:, :, i),  ce,  0.25);
end
conTime(4) = toc;

% Proposed B-spline.
tic
for i = 1 : numItr
    BsplineFittingFast( rawData(:, :, i), p, ce, er );
end
conTime(5) = toc;