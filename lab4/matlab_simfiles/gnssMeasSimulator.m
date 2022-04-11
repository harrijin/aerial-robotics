function [rpGtilde,rbGtilde] = gnssMeasSimulator(S,P)
% gnssMeasSimulator : Simulates GNSS measurements for quad.
%
%
% INPUTS
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
%
% OUTPUTS
%
% rpGtilde --- 3x1 GNSS-measured position of the quad's primary GNSS antenna,
%              in ECEF coordinates relative to the reference antenna, in
%              meters.
%
% rbGtilde --- 3x1 GNSS-measured position of secondary GNSS antenna, in ECEF
%              coordinates relative to the primary antenna, in meters.
%              rbGtilde is constrained to satisfy norm(rbGtilde) = b, where b
%              is the known baseline distance between the two antennas.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Harrison Jin
%+==============================================================================+  
% Constants
R_IG = Recef2enu(P.sensorParams.r0G);
rpB = P.sensorParams.raB(:,1);
rsB = P.sensorParams.raB(:,2);
% rpG
rpI = S.statek.rI+(S.statek.RBI)'*rpB;
rpG = (R_IG)'*rpI;
R_pG = R_IG'*P.sensorParams.RpL*R_IG;
omega_pG_cov = diag(diag(R_pG));
omega_pG = mvnrnd(zeros(3,1), omega_pG_cov)';
rpGtilde = rpG + omega_pG;
% rbG
rsI = S.statek.rI+(S.statek.RBI)'*rsB;
rsG = (R_IG)'*rsI;
rbG = rsG - rpG;
rbG_u = rbG/norm(rbG);
epsilon = 10^-8;
omega_bG_cov = norm(rbG)^2*P.sensorParams.sigmab^2*(eye(3,3)-(rbG_u*rbG_u'))+epsilon*eye(3,3);
omega_bG = mvnrnd(zeros(3,1), omega_bG_cov)';
rbGtilde = rbG + omega_bG;
rbGtilde = rbGtilde/norm(rbGtilde)*norm(rbG);
end