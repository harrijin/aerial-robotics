function [uCross] = crossProductEquivalent(u)
% crossProductEquivalent : Outputs the cross-product-equivalent matrix uCross
% such that for arbitrary 3-by-1 vectors u and v,
% cross(u,v) = uCross*v.
%
% INPUTS
%
% u ---------- 3-by-1 vector
%
%
% OUTPUTS
%
% uCross ----- 3-by-3 skew-symmetric cross-product equivalent matrix
%
%+------------------------------------------------------------------------------+
% References:
% Lecture Notes
%
% Author: Harrison Jin
%+==============================================================================+
    uCross = [     0  -1*u(3)    u(2); 
                 u(3)       0 -1*u(1);
              -1*u(2)     u(1)     0];
end

