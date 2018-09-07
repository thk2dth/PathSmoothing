function basisMatrix = MatrixBspline( curveParameter, degree, knotVector )
% obtain the matrix of the bspline basis function  with the given
% parameter.
% input:
%   parameter, the parameters where the basis matrix is derieved.
%                    could be a scalar or a vector.
%   degree, degree of the B-spline basis.
%   knotVector, knot vector of the B-spline basis.
% ouput:
%   basisMatrix, matrix of the basis function with a dimension of
%                     ( number of parameters * number of control points)

% this function wil call the function basisfun(i,u,p,U).
% HJ, 20160606
validateattributes(curveParameter, {'numeric'}, {'real', 'vector'});
validateattributes(degree, {'numeric'}, {'positive', 'integer', 'scalar'});
validateattributes(knotVector, {'numeric'}, {'real', 'vector', 'nondecreasing'});
numberKnotVector = length(knotVector);
numberControlPoints = numberKnotVector - degree - 1;
numberParameter = length(curveParameter);
basisMatrix = zeros(numberParameter, numberControlPoints);
for i = 1:numberParameter
    % in function findspan(), returned index is consistent with the that in "The NURBS Book". 
    index = findspan(numberControlPoints-1, degree, curveParameter(i), knotVector);
    tempN = basisfun(index, curveParameter(i), degree, knotVector);
    basisMatrix(i, index-degree+1 : index+1) = tempN;
end
end

