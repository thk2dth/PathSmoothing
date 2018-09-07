function [ctrlNew, knotsNew] = BsplineSplitting( degree, ctrl, knots, curveParameter )
% split the bspline by the curve parameters.
% after knot insert the multipilicity of each curve parameter will be p,
% where p is the degree of the bspline curve.
% input:
%   degree, degree of the bspline curve before splitting.
%   ctrl, control points of the bspline curve before inserting.
%   knots, knot vector of the bspline curve before inserting.
%   curveParameter, the curve parameter where the curve is splitted.
% output:
%   ctrlNew, the bspline after splitting.
%   knotsNew, knot vector after splitting.
% HJ, 20160607.
% Insert four intermediate knots between [u_i, u_{i+1}] rather than
% two knots.
% HJ, 20180116.
validateattributes(curveParameter, {'numeric'}, {'real', 'vector', 'nondecreasing'});
paraMul = FindMultiplicity(knots, curveParameter);
knotInsert = [];
lenPara = length(curveParameter);
for i = 2:lenPara-1 % the first and last curve parameters do not need to be inserted.
    knotInsert = [knotInsert, ones(1, degree-paraMul(i))*curveParameter(i)]; %#ok<*AGROW>
end
% firstly insert the curve parameter.
[ctrlNew, knotsNew] = bspkntins(degree, ctrl, knots, knotInsert);
% then insert the middle parameters if needed.
knotInsert = [];
indexStart = 2;
for i = 2:lenPara
    us = curveParameter(i-1);
    ue = curveParameter(i);
    indexEnd = find(abs(knotsNew-ue) < 1e-8, 1);
    num = indexEnd - indexStart;
    % insert four intermediate knots.
    if degree == num % N=0. Inserted four knots are evenly distributed between u_i and u_{i+1}.
        du = 0.2 * (ue - us);
        knotInsert = [knotInsert, us + du * [1, 2, 3, 4] ]; 
    elseif degree+1 == num % N=1. There are two situations.
        ul = knotsNew(indexEnd-1); % existing knot between u_i and u_{i+1}.
        um = 0.5 * (us + ue);
        if ul < um % situation 1.
            du = 0.25 * (ue - ul);
            knotInsert = [knotInsert, ul + du * [1, 2, 3] ];
        else % situation 2.
            du = 0.25 * (ul - us);
            knotInsert = [knotInsert, us + du * [1, 2, 3] ];
        end
    else % N=2. There are four situations.
        ul = knotsNew(indexEnd-2);
        ur = knotsNew(indexEnd-1);
        um = 0.5 * (us + ue);
        if ur <= um % situation 1.
            du = (ue - ur) / 3.0;
            knotInsert = [knotInsert, ur + du * [1, 2] ];
        elseif ul >= um % situation 2.
            du = (ul - us) / 3.0;
            knotInsert = [knotInsert, us + du * [1, 2] ];
        elseif ur - ul >= (ue - us) / 3.0 % situation 3.
            du = (ur - ul) / 3.0;
            knotInsert = [knotInsert, ul + du * [1, 2] ];
        else % situation 4.
            knotInsert = [knotInsert, 0.5 * (us + ul), 0.5 * (ur + ue)];
        end
    end
    indexStart = indexEnd;
end
[ctrlNew, knotsNew] = bspkntins(degree, ctrlNew, knotsNew, knotInsert);
end


function n = FindMultiplicity(knotVector, curveParameter)
% find the multiplicity of the curve parameter in knot vector.
% input:
%   knotVector, knot vector.
%   curveParameter, the parameter to find the multiplicity
% output:
%   n, the multiplicity of curve parameter in knot vector.
% HJ, 20160607.
validateattributes(knotVector, {'numeric'}, {'real', 'vector', 'nondecreasing'});
validateattributes(curveParameter, {'numeric'}, {'real', 'vector', 'nondecreasing'});
lenPara = length(curveParameter);
n = zeros(1, lenPara);
lenKnot = length(knotVector);
for i = 1:lenPara
    for j = 1:lenKnot
        if abs(knotVector(j) - curveParameter(i)) < 1e-6
            n(i) = n(i) + 1;
        end
    end
end
end