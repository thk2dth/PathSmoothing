function d = Distance2Line( startPoint, endPoint, evaluatePoint )
% distance from evaluatePoint to the line which crosses startPoint and
% endPoint.
% calculate the height of the triangle to obtain the distance
% input:
%   startPoint, endPoint,  two points on the line.
%   evaluatePoint, the point of which the distance to be evaluated
% output:
%   d, the distance from evaluatePoint to the line
% HJ, 20160607.
validateattributes(startPoint, {'numeric'}, {'real', 'column'});
validateattributes(endPoint, {'numeric'}, {'real', 'column'});
validateattributes(evaluatePoint, {'numeric'}, {'real', 'column'});
v1 = startPoint - evaluatePoint;
v2 = endPoint - evaluatePoint;
if 2 == length(v1)
    v1 = [v1; 0];
    v2 = [v2; 0];
end
v3 = endPoint - startPoint;
d = norm(cross(v1, v2)) / norm(v3);
end