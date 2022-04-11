function [Q] = simulateQuadrotorControl(R,S,P)
% simulateQuadrotorControl : Simulates closed-loop control of a quadrotor
%                            aircraft.
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%          tVec = Nx1 vector of uniformly-sampled time offsets from the
%                 initial time, in seconds, with tVec(1) = 0.
%
%        rIstar = Nx3 matrix of desired CM positions in the I frame, in
%                 meters.  rIstar(k,:)' is the 3x1 position at time tk =
%                 tVec(k).
%
%        vIstar = Nx3 matrix of desired CM velocities with respect to the I
%                 frame and expressed in the I frame, in meters/sec.
%                 vIstar(k,:)' is the 3x1    velocity at time tk = tVec(k).
%
%        aIstar = Nx3 matrix of desired CM accelerations with respect to the I
%                 frame and expressed in the I frame, in meters/sec^2.
%                 aIstar(k,:)' is the 3x1 acceleration at time tk =
%                 tVec(k).
%
%        xIstar = Nx3 matrix of desired body x-axis direction, expressed as a
%                 unit vector in the I frame. xIstar(k,:)' is the 3x1
%                 direction at time tk = tVec(k).
%  
% S ---------- Structure with the following elements:
%
%  oversampFact = Oversampling factor. Let dtIn = R.tVec(2) - R.tVec(1). Then
%                 the output sample interval will be dtOut =
%                 dtIn/oversampFact. Must satisfy oversampFact >= 1.
%
%        state0 = State of the quad at R.tVec(1) = 0, expressed as a structure
%                 with the following elements:
%                   
%                   r = 3x1 position in the world frame, in meters
% 
%                   e = 3x1 vector of Euler angles, in radians, indicating the
%                       attitude
%
%                   v = 3x1 velocity with respect to the world frame and
%                       expressed in the world frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
%       distMat = (N-1)x3 matrix of disturbance forces acting on the quad's
%                 center of mass, expressed in Newtons in the world frame.
%                 distMat(k,:)' is the constant (zero-order-hold) 3x1
%                 disturbance vector acting on the quad from R.tVec(k) to
%                 R.tVec(k+1).
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%  sensorParams = Structure containing sensor parameters, as defined in
%                 sensorParamsScript.m
%
%
% OUTPUTS
%
% Q ---------- Structure with the following elements:
%
%          tVec = Mx1 vector of output sample time points, in seconds, where
%                 Q.tVec(1) = R.tVec(1), Q.tVec(M) = R.tVec(N), and M =
%                 (N-1)*oversampFact + 1.
%  
%         state = State of the quad at times in tVec, expressed as a
%                 structure with the following elements:
%                   
%                rMat = Mx3 matrix composed such that rMat(k,:)' is the 3x1
%                       position at tVec(k) in the I frame, in meters.
% 
%                eMat = Mx3 matrix composed such that eMat(k,:)' is the 3x1
%                       vector of Euler angles at tVec(k), in radians,
%                       indicating the attitude.
%
%                vMat = Mx3 matrix composed such that vMat(k,:)' is the 3x1
%                       velocity at tVec(k) with respect to the I frame
%                       and expressed in the I frame, in meters per
%                       second.
%                 
%           omegaBMat = Mx3 matrix composed such that omegaBMat(k,:)' is the
%                       3x1 angular rate vector expressed in the body frame in
%                       radians, that applies at tVec(k).
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  
% Load parameters
params.quadParams = P.quadParams;
params.constants = P.constants;
% Number of iterations to run
N = length(R.tVec);
% Sampling rate
dtOut = (R.tVec(2)-R.tVec(1))/S.oversampFact;
% Create empty storage vectors
tVecOut = [];
rMat = [];
vMat = [];
eMat = [];
omegaBMat = [];

% Set initial states
rK = S.state0.r;
vK = S.state0.v;
eK = S.state0.e;
omegaB_K = S.state0.omegaB;
% Initial hover state
omegaK = [0 0 0 0]';

% Iterate
for k = 1:N-1
    tspan = [R.tVec(k):dtOut:R.tVec(k+1)]';
    % Build X_big
    RBI = euler2dcm(eK);
    Xk = [rK;vK;RBI(1:3,1);RBI(1:3,2);RBI(1:3,3);omegaB_K;omegaK];
    % Initialize ode45 inputs
    distVec = S.distMat(k,:)';
    % Find Fk and zIstark using trajectory controller
    Rk.rIstark = R.rIstar(k,:)';
    Rk.vIstark = R.vIstar(k,:)';
    Rk.aIstark = R.aIstar(k,:)';
    Sk.statek.rI = rK;
    Sk.statek.RBI = RBI;
    Sk.statek.vI = vK;
    Sk.statek.omegaB = omegaB_K;
    [Fk, Rk.zIstark] = trajectoryController(Rk,Sk,P);
    % Find Nbk using attitude controller
    Rk.xIstark = R.xIstar(k,:)';
    NBk = attitudeController(Rk,Sk,P);
    % Find desired voltages
    eaVec = voltageConverter(Fk,NBk,P);
    % Run ODE solver on segment
    [tVecK,XMatk] = ode45(@(t,X)quadOdeFunctionHF(t,X,eaVec, distVec, params), tspan, Xk);
    % Store outputs
    tVecOut = [tVecOut; tVecK(1:end-1)];
    rMat = [rMat; XMatk(1:end-1, 1:3)];
    vMat = [vMat; XMatk(1:end-1, 4:6)];
    for row = 1:size(XMatk,1)-1
        rotMat = [XMatk(row,7:9); XMatk(row,10:12); XMatk(row,13:15)]';
        eMat = [eMat;(dcm2euler(rotMat))'];
    end
    omegaBMat = [omegaBMat;XMatk(1:end-1, 16:18)];
    % Prep for next iteration
    rK = XMatk(end, 1:3)';
    vK = XMatk(end, 4:6)';
    eK = dcm2euler([XMatk(end,7:9); XMatk(end,10:12); XMatk(end,13:15)]');
    omegaB_K = XMatk(end, 16:18)';
    omegaK = XMatk(end, 19:22)';
end
% Store the final state of the final segment
tVecOut = [tVecOut; tVecK(end,:)];
rMat = [rMat; rK'];
vMat = [vMat; vK'];
eMat = [eMat; eK'];
omegaBMat = [omegaBMat; omegaB_K'];

% Store into output structure
state.rMat = rMat;
state.vMat = vMat;
state.eMat = eMat;
state.omegaBMat = omegaBMat;
Q.tVec = tVecOut;
Q.state = state;
end