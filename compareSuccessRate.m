% Compare the success rate of the different fitting methods.
% The comparison consists of the following steps.
% (1) Generate 5*100 points in 3-D space. The distance and angle between the adjacent
% points are both randomly generated and belong to a specified range.
% (2) The Akima fittting, the fast Bezier fitting, the iterative B-spline
% fitting, and the proposed two-step methods are applied to the generated
% 500 points. Hence there is 100 loops, and each loop processes 5 points.
% (3) 1000 points are sampled on the resultant curves, and are used to
% calculate the chord errors.
% HJ, 20180123.

clear all;
numItr = 100;
numPts = 5;
ce = 0.04; % chord error
er = 0.2; % error ratio. The smaller, the more strigent for rough fitting.
p = 3; % degree of B-spline.
lenRange = [0.5, 2];
alphaRange = [2.5, 30] * (pi/180); % the angle range between adjacent lines.
betaRange = [-pi, pi]; % no limit for beta.
numInp = 1000;
numBiChordErrorPass = 0; % Data points that pass the bi-chord error test.
numFail = zeros(1, 5);
numIterativeBspline = 50;
rawData = zeros(3, numPts, numItr);
maxErr = zeros(numItr, 5);
wl = 2; % wide line

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
    tf = 0;
    % Bi-chord error test.
    for k = 2 : numPts-1
        [flag, ~, ~] = BiChordErrorTest(rawData(:, k-1, i), rawData(:, k, i),...
            rawData(:, k+1, i), ce);
        tf = tf + flag;
    end
    if tf == numPts - 2 % All points pass the test.
        numBiChordErrorPass = numBiChordErrorPass + 1;
    end
end

u = linspace(0, 1, numInp); % sampled curve parameters.
for i = 1 : numItr
    % Akima fitting.
    [ctrl, knots] = AkimaFitting5Points( rawData(:, :, i) );
    for j = 1 : numInp
        pt = AkimaPoints(knots, ctrl, u(j) );
        conErr = contourErrorPolygonDirect(pt', rawData(:, :, i)');
        if conErr > maxErr(i, 1)
            maxErr(i, 1) = conErr;
        end
    end
    if maxErr(i, 1) > ce
        numFail(1) = numFail(1) + 1;
    end
    
    % Fast Bezier fitting.
    [ctrl, knots] = BezierFittingFast( rawData(:, :, i) );
    pts = BezierFastPoints(knots, ctrl, numInp);
    for j = 1 : numInp
        conErr = contourErrorPolygonDirect(pts(:, j)', rawData(:, :, i)');
        if conErr > maxErr(i, 2)
            maxErr(i, 2) = conErr;
        end
    end
    if maxErr(i, 2) > ce
        numFail(2) = numFail(2) + 1;
    end
    
    % Iterative B-spline fitting.
    [ctrl, knots] = BsplineFittingIterative( rawData(:, :, i), p, ce, numIterativeBspline);
    for j = 1 : numInp
        pt = BsplinePoints(p, ctrl, knots, u(j) );
        conErr = contourErrorPolygonDirect(pt', rawData(:, :, i)');
        if conErr > maxErr(i, 3)
            maxErr(i, 3) = conErr;
        end
    end
    if maxErr(i, 3) > ce
        numFail(3) = numFail(3) + 1;
    end
    
    % Transition
    [ctrl, knots, p] = BsplineTransition(rawData(:, :, i), ce, 0.25);
    for j = 1 : numInp
        pt = BsplinePoints(p, ctrl, knots, u(j) );
        conErr = contourErrorPolygonDirect(pt', rawData(:, :, i)');
        if conErr > maxErr(i, 4)
            maxErr(i, 4) = conErr;
        end
    end
    if maxErr(i, 4) > ce
        numFail(4) = numFail(4) + 1;
    end
    
    % Proposed B-spline fitting.
    [ctrl, knots] = BsplineFittingFast( rawData(:, :, i), p, ce, er );
    for j = 1 : numInp
        pt = BsplinePoints(p, ctrl, knots, u(j) );
        conErr = contourErrorPolygonDirect(pt', rawData(:, :, i)');
        if conErr > maxErr(i, 5)
            maxErr(i, 5) = conErr;
        end
    end
    if maxErr(i, 5) > ce
        numFail(5) = numFail(5) + 1;
    end
    
end

successRate = 1 -  numFail / numItr;

%% Figure
figure('Name', 'Maximum error')
subplot(3, 2, 1)
xAxis = 1 : numItr;
plot(xAxis, maxErr(:, 1), 'ro');
hold on;
plot([1, numItr], [ce, ce], 'k-.', 'LineWidth', wl);
% xlabel('{Instance}');
% ylabel('{\bfError} (mm)');

subplot(3, 2, 2)
plot(xAxis, maxErr(:, 2), 'gd');
hold on;
plot([1, numItr], [ce, ce], 'k-.', 'LineWidth', wl);
% xlabel('{Instance}');
% ylabel('{\bfError} (mm)');

subplot(3, 2, 3)
plot(xAxis, maxErr(:, 3), 'b+');
hold on;
plot([1, numItr], [ce, ce], 'k-.', 'LineWidth', wl);
% xlabel('{Instance}');
% ylabel('{\bfError} (mm)');

subplot(3, 2, 4)
plot(xAxis, maxErr(:, 4), 'c*');
hold on;
plot([1, numItr], [ce, ce], 'k-.', 'LineWidth', wl);
% xlabel('{Instance}');
% ylabel('{\bfError} (mm)');

subplot(3, 2, [5, 6])
plot(xAxis, maxErr(:, 5), 'ms');
hold on;
plot([1, numItr], [ce, ce], 'k-.', 'LineWidth', wl);
% xlabel('{Instance}');
% ylabel('{\bfError} (mm)');
