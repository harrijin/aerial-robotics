function [NBk] = attitudeController(R,S,P)
% attitudeController : Controls quadcopter toward a reference attitude
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%       zIstark = 3x1 desired body z-axis direction at time tk, expressed as a
%                 unit vector in the I frame.
%
%       xIstark = 3x1 desired body x-axis direction, expressed as a
%                 unit vector in the I frame.
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude
%
%                  vI = 3x1 velocity with respect to the I frame and
%                       expressed in the I frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%
% OUTPUTS
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Harrison Jin
%+==============================================================================+  
% Initialize controller gains
K = diag([1 1.5 0.1]);
KD = diag([0.25 0.3 0.05]);

% Find desired attitude
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];
yIstark = crossProductEquivalent(R.zIstark)*R.xIstark;
yIstark = yIstark/norm(yIstark);
zIstark = R.zIstark;
xIstark = crossProductEquivalent(yIstark)*zIstark;
RBIstark = [xIstark, yIstark, zIstark]';
% Find error matrix and eigenvector
rE = RBIstark*(S.statek.RBI)';
eE = [rE(2,3) - rE(3,2);
      rE(3,1) - rE(1,3);
      rE(1,2) - rE(2,1)];
% Find needed torque
NBk = K*eE - KD*S.statek.omegaB + crossProductEquivalent(S.statek.omegaB)*P.quadParams.Jq*S.statek.omegaB;
end