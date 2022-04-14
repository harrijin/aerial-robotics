function [rXIHat,Px] = estimate3dFeatureLocation(M,P)
% estimate3dFeatureLocation : Estimate the 3D coordinates of a feature point
%                             seen by two or more cameras with known pose.
%
%
% INPUTS
%
% M ---------- Structure with the following elements:
%
%       rxArray = 1xN cell array of measured positions of the feature point
%                 projection on the camera's image plane, in pixels.
%                 rxArray{i} is the 2x1 vector of coordinates of the feature
%                 point as measured by the ith camera.  To ensure the
%                 estimation problem is observable, N must satisfy N >= 2 and
%                 at least two cameras must be non-colinear.
%
%      RCIArray = 1xN cell array of I-to-camera-frame attitude matrices.
%                 RCIArray{i} is the 3x3 attitude matrix corresponding to the
%                 measurement rxArray{i}.
%
%       rcArray = 1xN cell array of camera center positions.  rcArray{i} is
%                 the 3x1 position of the camera center corresponding to the
%                 measurement rxArray{i}, expressed in the I frame in meters.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the quad's
%                 sensors, as defined in sensorParamsScript.m
%
% OUTPUTS
%
%
% rXIHat -------- 3x1 estimated location of the feature point expressed in I
%                 in meters.
%
% Px ------------ 3x3 error covariance matrix for the estimate rxIHat.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Harrison Jin
%+==============================================================================+ 
    H=[];
    for i = 1:size(M.rxArray,2)
        proj_mat = P.sensorParams.K*[M.RCIArray{i}, -1*M.RCIArray{i}*M.rcArray{i}];
        x_tilde = M.rxArray{i}(1)*P.sensorParams.pixelSize;
        y_tilde = M.rxArray{i}(2)*P.sensorParams.pixelSize;
        H=[H;
           x_tilde*proj_mat(3,:)-proj_mat(1,:);
           y_tilde*proj_mat(3,:)-proj_mat(2,:)];
    end
    Hr = H(:,1:3);
    z = -1*H(:,4);
    RcArray = repmat({P.sensorParams.Rc}, 1, size(M.rxArray,2));
    R=P.sensorParams.pixelSize^2*blkdiag(RcArray{:});
    Px=inv(Hr'*inv(R)*Hr);
    rXIHat=Px*Hr'*inv(R)*z;
end