function [e] = dcm2euler(R_BW)
% dcm2euler : Converts a direction cosine matrix R_BW to Euler angles phi =
%             e(1), theta = e(2), and psi = e(3) (in radians) for a 3-1-2
%             rotation. If the conversion to Euler angles is singular (not
%             unique), then this function issues an error instead of returning
%             e.
%
% Let the world (W) and body (B) reference frames be initially aligned.  In a
% 3-1-2 order, rotate B away from W by angles psi (yaw, about the body Z
% axis), phi (roll, about the body X axis), and theta (pitch, about the body Y
% axis).  R_BW can then be used to cast a vector expressed in W coordinates as
% a vector in B coordinates: vB = R_BW * vW
%
% INPUTS
%
% R_BW ------- 3-by-3 direction cosine matrix 
%
%
% OUTPUTS
%
% e ---------- 3-by-1 vector containing the Euler angles in radians: phi =
%              e(1), theta = e(2), and psi = e(3).  By convention, these
%              should be constrained to the following ranges: -pi/2 <= phi <=
%              pi/2, -pi <= theta < pi, -pi <= psi < pi.  
% 
%+------------------------------------------------------------------------------+
% References:
% Lecture notes
%
% Author: Harrison Jin
%+==============================================================================+  
% Singular case: phi == pi/2 or multiple
    if abs(R_BW(2,3) - 1) < 0.001
        error('Input matrix is singular case, unable to convert to euler angles');
    end
% Nominal case
    phi = asin(R_BW(2,3));
    psi = atan2(-1*R_BW(2,1), R_BW(2,2));
    theta = atan2(-1*R_BW(1,3), R_BW(3,3));
    e = [phi;theta;psi];
end