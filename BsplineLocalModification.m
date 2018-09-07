function [ w, ctrl ] = BsplineLocalModification( us, ue, pts, pte, der1s, der2s, der1e, der2e, me )
% locally modify the bspline under the given constraints.
% a optimization model which makes the inserted knots
% evenly distributed in the interval [u1, u2) is established.
% input:
%   us, ue, start/end curve parameter.
%   pts, pte, start/end control points
%   der1s, der1e, first derivatives at the start/end parameter
%   der2s, der2e, second derivatives at the start/end parameter.
%   me, modification error.
% output:
%   w, knots used to replace the original ones.
%   ctrl, control points used to replace the original ones.
% HJ, 20160608.
% us < w1 < w2 < w3 < w4 < ue
validateattributes(us , {'numeric'}, {'real', 'scalar'});
validateattributes(ue , {'numeric'}, {'real', 'scalar'});
validateattributes(pts, {'numeric'}, {'real', 'column'});
validateattributes(pte, {'numeric'}, {'real', 'column'});
validateattributes(der1s, {'numeric'}, {'real', 'column'});
validateattributes(der2s, {'numeric'}, {'real', 'column'});
validateattributes(der1e, {'numeric'}, {'real', 'column'});
validateattributes(der2e, {'numeric'}, {'real', 'column'});
validateattributes(me, {'numeric'}, {'real', 'scalar'});

boundaryCondition = [pts, pte, der1s, der1e, der2s, der2e];
w = CalculateWeight(us, ue, boundaryCondition, me);
S8 = pts + der1s*(w(1)-us) / 3.0;
S13 = pte - der1e*(ue-w(4)) / 3.0;
S9 = pts + der1s*(w(1)+w(2)-2*us) / 3.0 + der2s*(w(1)-us)*(w(2)-us) / 6.0;
S12 = pte - der1e*(2*ue-w(3)-w(4)) / 3.0 + der2e*(ue-w(3))*(ue-w(4)) / 6.0;
S10 = (2*S9 + S12) / 3.0;
S11 = (S9 + 2*S12) / 3.0;
ctrl = [S8, S9, S10, S11, S12, S13];
end


function w = CalculateWeight( us, ue, boundaryCondition, me )
% weight calculation by given chordError.
% input:
%   us, ue, start/end curve parameter.
%   boundaryCondition, start/end points, start/end first order derivatives,
%                                  start/end second order derivatives.
%   em, fitting error.
% output:
%   w, w0, w1, w2, w3.
t = boundaryCondition(:, 2) - boundaryCondition(:, 1); % line R7R14
n12 = cross(boundaryCondition(:, 3), t);
n12 = n12 / norm( n12 );
n11 = cross(n12, t);
n11 = n11 / norm( n11 );
n22 = cross(boundaryCondition(:, 4), t);
n22 = n22 / norm(n22);
n21 = cross(n22, t);
n21 = n21 / norm(n21);
c1 = dot(n11, boundaryCondition(:, 3) );
c2 = dot(n12, boundaryCondition(:, 3) ); % c2 = 0
c3 = dot(n11, boundaryCondition(:, 5) );
c4 = dot(n12, boundaryCondition(:, 5) );
c5 = dot(n21, boundaryCondition(:, 4) );
c6 = dot(n22, boundaryCondition(:, 4) ); % c6 = 0
c7 = dot(n21, boundaryCondition(:, 6) );
c8 = dot(n22, boundaryCondition(:, 6) );
du = 0.2 * (ue - us);
[dw1, dw2] = CalculateLambda(c1, c3, c4, du, me);
w(1) = us + dw1;
w(2) = us + dw2;
[dw4, dw3] = CalculateLambda(-c5, c7, c8, du, me);
w(4) = ue - dw4;
w(3) = ue - dw3;
end

function [lamb0, lamb1] = CalculateLambda(c1, c3, c4, du, em)
if c3 == 0
    x1 = 1.5 * em / abs( c1 );
elseif c3 < 0
    x1 = du;
else % c3 > 0
    if c1 > 0 % c1 > 0
        x1 = (-2*c1 + sqrt(4*c1*c1+6*c3*em) ) / c3;
    else
        if c1*c1 > 1.5*c3*em % c1 > 0
            x1 = (-2*c1 - sqrt(4*c1*c1-6*c3*em) ) / c3;
        else % 2*c1*c1<=3*c3*em, c1 <= 0
            x1 = (-2*c1 + sqrt(4*c1*c1+6*c3*em) ) / c3;
        end
    end
end
lamb0 = min(min(du, 1.5*em/abs(c1) ), x1);
k1 = c1 * lamb0 / 3;
k2 = c1 / 3 + c3 * lamb0 / 6;
k3 = c4 * lamb0 / 6;
x2 = (-k1*k2 + sqrt(k2*k2*em*em + k3*k3*(em*em - k1*k1) ) ) / (k2*k2 + k3*k3);
lamb1 = min(2*du, x2);
end