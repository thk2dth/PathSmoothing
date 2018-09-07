% Compare different fitting methods.
% The Akima fitting and fast Bezier fitting methods do
% not constrain the chord error directly.
len = 1;
% rawData =  [0, 1, 2, 3, 4; 0, 0.2, 0, -0.2, 0; 0, 0, 0, 0, 0] * 0.5; % eg 1.
rawData =  [0, 1, 2, 3, 4; 0, 0.1, 0, -0.12, 0; 0, 0, 0, 0, 0]; % eg 2.
N = 2000; % number of points to be evaluated.
ce = 0.04; % chord error.
wl = 2;


%% Akima curve fitting
[controlPoints, knotVector ] = AkimaFitting5Points( rawData );
u = linspace(0, 1, N);
[pts, ~, ~, cur]  = AkimaPoints(knotVector, controlPoints, u);
err = contourErrorPolygonDirect(pts', rawData');
fc = figure('Name', 'Contour after fitting');
plot(pts(1,: ), pts(2, :), 'r' );
xlabel('{\bfX}(mm)');
ylabel('{\bfY}(mm)');
hold on;

fr = figure('Name', 'Contour error');
plot(pts(1, :), err, 'r'); 
xlabel('{\bfX}(mm)');
ylabel('{\bfChord error}(mm)');
hold on;

fcur = figure('Name', 'Curvature');
plot(pts(1, :), cur, 'r'); 
xlabel('{\bfX}(mm)');
ylabel('{\bfCurvature}(mm^{-1})');
hold on;

%% Fast Bezier curve fitting
[controlPoints, knots] = BezierFittingFast( rawData );
[pts, ~, ~, cur] = BezierFastPoints(knots, controlPoints, N);
figure(fc);
plot(pts(1, :), pts(2, :), 'g-' );

err = contourErrorPolygonDirect(pts', rawData');
figure(fr)
plot(pts(1, :), err, 'g'); 

figure(fcur)
plot(pts(1, :), cur, 'g'); 

%% Iterative B-spline fitting
p = 3; % degree of the B-spline.
[ctrl, knots] = BsplineFittingIterative( rawData, p, ce, 50 );
u = linspace(0, 1, N);
[pts, ~, ~, cur]  = BsplinePoints(p, ctrl, knots, u );
figure(fc);
plot(pts(1, :), pts(2, :), 'b-' );

err = contourErrorPolygonDirect(pts', rawData');
figure(fr)
plot(pts(1, :), err, 'b'); 

figure(fcur)
plot(pts(1, :), cur, 'b'); 

%% Transition
[ctrl, knots, p] = BsplineTransition(rawData, ce, 0.25);
u = linspace(0, 1, N);
[pts, ~, ~, cur]  = BsplinePoints(p, ctrl, knots, u );
err = contourErrorPolygonDirect(pts', rawData');

figure(fc)
plot(pts(1,: ), pts(2, :), 'c' );

figure(fr)
plot(pts(1, :), err, 'c'); 

figure(fcur)
plot(pts(1, :), cur, 'c'); 




%% Proposed method
er = 0.2;
degree = 3;
[ctrl, knots] = BsplineFittingFast( rawData, degree, ce, er );
u = linspace(0, 1, N);
[pts, der1, der2, cur]  = BsplinePoints(degree, ctrl, knots, u);
figure(fc);
plot(pts(1, :), pts(2, :), 'm-', 'LineWidth', wl );
figure(fc);
plot(rawData(1, :), rawData(2, :), 'ko-.', 'LineWidth', wl);
legend('Akima', 'Bezier', 'Iterative', 'Transition',...
    'Proposed', 'Linear');
% xlim([-0.02, 2]);

err = contourErrorPolygonDirect(pts', rawData');
figure(fr)
plot(pts(1, :), err, 'm', 'LineWidth', wl); 
plot([pts(1, 1), pts(1, end)], [ce, ce], 'k-.', 'LineWidth', wl);
legend('Akima', 'Bezier', 'Iterative', 'Transition', ...
    'Proposed', 'Tolerance');
% xlim([0, 2]);

figure(fcur)
plot(pts(1, :), cur, 'm-', 'LineWidth', wl ); 
legend('Akima', 'Bezier', 'Iterative', 'Transition', ...
    'Proposed');
% xlim([0, 2]);

% derivatives.
der = zeros(2, N);
for i = 1 : N
    der(1, i) = norm(der1(:, i) );
    der(2, i) = norm(der2(:, i) );
end
fd = figure('Name', 'Derivative');
plotyy(pts(1, :), der(1, :), pts(1, :), der(2, :) );