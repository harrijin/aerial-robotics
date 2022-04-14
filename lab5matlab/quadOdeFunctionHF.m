function [Xdot] = quadOdeFunctionHF(t,X,eaVec,distVec,P)
% quadOdeFunctionHF : Ordinary differential equation function that models
%                     quadrotor dynamics -- high-fidelity version.  For use
%                     with one of Matlab's ODE solvers (e.g., ode45).
%
%
% INPUTS
%
% t ---------- Scalar time input, as required by Matlab's ODE function
%              format.
%
% X ---------- Nx-by-1 quad state, arranged as 
%
%              X = [rI',vI',RBI(1,1),RBI(2,1),...,RBI(2,3),RBI(3,3),...
%                   omegaB',omegaVec']'
%
%              rI = 3x1 position vector in I in meters
%              vI = 3x1 velocity vector wrt I and in I, in meters/sec
%             RBI = 3x3 attitude matrix from I to B frame
%          omegaB = 3x1 angular rate vector of body wrt I, expressed in B
%                   in rad/sec
%        omegaVec = 4x1 vector of rotor angular rates, in rad/sec.
%                   omegaVec(i) is the angular rate of the ith rotor.
%
%    eaVec --- 4x1 vector of voltages applied to motors, in volts.  eaVec(i)
%              is the constant voltage setpoint for the ith rotor.
%
%  distVec --- 3x1 vector of constant disturbance forces acting on the quad's
%              center of mass, expressed in Newtons in I.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
% OUTPUTS
%
% Xdot ------- Nx-by-1 time derivative of the input vector X
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  Harrison Jin
%+==============================================================================+
% Extract data from input X
vI = X(4:6);
R_BI = zeros(3,3);
R_BI(1:3,1) = X(7:9);
R_BI(1:3,2) = X(10:12);
R_BI(1:3,3) = X(13:15);
omegaB = X(16:18);
omegaVec = X(19:22);
% Get forces and torques from rotors
torqueDirs = -1*P.quadParams.kN' .* P.quadParams.omegaRdir;
rotorForces = P.quadParams.kF' *(omegaVec.^2);
rotorTorques = torqueDirs *(omegaVec.^2);
% Initialize Xdot
Xdot = zeros(22,1);
% rI dot
Xdot(1:3) = vI;
% vI dot
weight = [0;0;P.quadParams.m * -1 * P.constants.g];
rotorThrust = R_BI' * [0;0;rotorForces];
zI = R_BI' * [0;0;1];
drag = 0;
if norm(vI) > 0
    vIu = vI/norm(vI);
    drag = -0.5*P.quadParams.Cd*P.quadParams.Ad*P.constants.rho*(zI'*vI*norm(vI))*vIu;
end
Xdot(4:6) = (weight + rotorThrust + distVec + drag)/P.quadParams.m;
% RBI dot
RBI_dot = -1*crossProductEquivalent(omegaB)*R_BI;
Xdot(7:9) = RBI_dot(1:3,1);
Xdot(10:12) = RBI_dot(1:3,2);
Xdot(13:15) = RBI_dot(1:3,3);
% omegaB dot
torque = [0;0;rotorTorques];
for k = 1:4
    % Torque due to thrust
    torque = torque + crossProductEquivalent(P.quadParams.rotor_loc(1:3,k)) ...
             * (P.quadParams.kF(k)*[0;0;omegaVec(k)^2]);
end
Xdot(16:18) = inv(P.quadParams.Jq)*(torque-crossProductEquivalent(omegaB)*P.quadParams.Jq*omegaB);
% omegaVec dot
Xdot(19:22) = (P.quadParams.cm.*eaVec - omegaVec)./(P.quadParams.taum);
end
