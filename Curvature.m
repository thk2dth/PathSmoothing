function cur = Curvature( d1, d2 )
% Calculate the absolute curvature according to the first
% and second order derivatives.
cur = norm(cross(d1, d2) ) / norm(d1)^3;
end

