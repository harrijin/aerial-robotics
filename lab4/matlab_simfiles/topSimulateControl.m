% Top-level script for calling simulateQuadrotorControl or
% simulateQuadrotorEstimationAndControl
% 'clear all' is needed to clear out persistent variables from run to run
clear all; clc;
%% Read trajectory data in
pos = readmatrix("trajectory/pos.csv")';
vel = readmatrix("trajectory/vel.csv")';
acc = readmatrix("trajectory/acc.csv")';
N = size(pos,1);
%%
% Seed Matlab's random number: this allows you to simulate with the same noise
% every time (by setting a nonnegative integer seed as argument to rng) or
% simulate with a different noise realization every time (by setting
% 'shuffle' as argument to rng).
rng('shuffle');
% Assert this flag to call the full estimation and control simulator;
% otherwise, only the control simulator is called
estimationFlag = 1;
% Populate reference trajectory
R.tVec = pos(:,1);
R.rIstar = [pos(:,2),pos(:,3),zeros(N,1)];
R.vIstar = [vel(:,2),vel(:,3),zeros(N,1)];
R.aIstar = [acc(:,2),acc(:,3),zeros(N,1)];
% Change trajectory to match coordinate system of occupancy grid
R.rIstar(:,2) = -1*R.rIstar(:,2);
R.vIstar(:,2) = -1*R.vIstar(:,2);
R.aIstar(:,2) = -1*R.aIstar(:,2);
% The desired xI points toward the origin. The code below also normalizes
% each row in R.xIstar.
R.xIstar = ([1 0 0]'*ones(1,N))';
% Matrix of disturbance forces acting on the body, in Newtons, expressed in I
S.distMat = 0*randn(N-1,3);
% Initial position in m
S.state0.r = [0 0 0]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e = [0 0 0]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v = [0 0 0]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
S.state0.omegaB = [0 0 0]';
% Oversampling factor
S.oversampFact = 2;
% Feature locations in the I frame
S.rXIMat = [0,0,0; 0,0,0.7]; 
S.rXIMat = [];
% Quadrotor parameters and constants
quadParamsScript;
constantsScript;
sensorParamsScript;
P.quadParams = quadParams; 
P.constants = constants; 
P.sensorParams = sensorParams;

if(estimationFlag)
  Q = simulateQuadrotorEstimationAndControl(R,S,P);
else
  Q = simulateQuadrotorControl(R,S,P);
end

S2.tVec = Q.tVec;
S2.rMat = Q.state.rMat;
S2.eMat = Q.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = false;
S2.gifFileName = 'testGif.gif';
S2.bounds=[0 10 -10 0 -1 1];
visualizeQuad(S2);

figure(2);clf;
plot(Q.tVec,Q.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 

figure(5);clf;
plot(Q.state.rMat(:,1), Q.state.rMat(:,2)); 
axis equal; grid on; hold on;
xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');



