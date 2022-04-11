function [rx] = hdCameraSimulator(rXI,S,P)
% hdCameraSimulator : Simulates feature location measurements from the
%                     quad's high-definition camera. 
%
%
% INPUTS
%
% rXI -------- 3x1 location of a feature point expressed in I in meters.
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position of CM in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude of B frame wrt I frame
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the
%                 quad's sensors, as defined in sensorParamsScript.m 
%
% OUTPUTS
%
% rx --------- 2x1 measured position of the feature point projection on the
%              camera's image plane, in pixels.  If the feature point is
%              not visible to the camera (the ray from the feature to the
%              camera center never intersects the image plane), then rx is
%              an empty matrix.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Harrison Jin
%+==============================================================================+  
RCI = P.sensorParams.RCB*S.statek.RBI;
t = -1*RCI*(S.statek.rI+S.statek.RBI'*P.sensorParams.rocB);
xc=[RCI t]*[rXI;1];
if xc(3) < 0
    rx = [];
else
    xc=P.sensorParams.K*xc;
    xc = xc/xc(3);
    if (abs(xc(1)) > P.sensorParams.imagePlaneSize(1)/2) || (abs(xc(2)) > P.sensorParams.imagePlaneSize(2)/2)
        rx = [];
    else
        xc = xc(1:2);
        omega_c = mvnrnd(zeros(2,1),P.sensorParams.Rc)';
        rx=(1/(P.sensorParams.pixelSize))*xc+omega_c;
    end
end
end