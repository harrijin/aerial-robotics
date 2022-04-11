% Top-level script for calling simulateQuadrotorControl or
% simulateQuadrotorEstimationAndControl

% 'clear all' is needed to clear out persistent variables from run to run
clear all; clc;
% Seed Matlab's random number: this allows you to simulate with the same noise
% every time (by setting a nonnegative integer seed as argument to rng) or
% simulate with a different noise realization every time (by setting
% 'shuffle' as argument to rng).
rng('shuffle');
% Assert this flag to call the full estimation and control simulator;
% otherwise, only the control simulator is called
estimationFlag = 1;
% Total simulation time, in seconds
Tsim = 10;
% Update interval, in seconds
delt = 0.005;
% Time vector, in seconds 
N = floor(Tsim/delt);
tVec=[0:N-1]'*delt;
% Angular rate of orbit, in rad/sec
n = 2*pi/Tsim;
% Radius of circle, in meters
r = 4;
% Populate reference trajectory
R.tVec = tVec;
R.rIstar = [r*cos(n*tVec),r*sin(n*tVec),ones(N,1)];
R.vIstar = [-r*n*sin(n*tVec),r*n*cos(n*tVec),zeros(N,1)];
R.aIstar = [-r*n*n*cos(n*tVec),-r*n*n*sin(n*tVec),zeros(N,1)];
% The desired xI points toward the origin. The code below also normalizes
% each row in R.xIstar.
R.xIstar = diag(1./vecnorm(R.rIstar'))*(-R.rIstar);
% Matrix of disturbance forces acting on the body, in Newtons, expressed in I
S.distMat = 0*randn(N-1,3);
% Initial position in m
S.state0.r = [r 0 0]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e = [0 0 pi]';
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
S2.bounds=4.5*[-1 1 -1 1 -0.1 1];
visualizeQuad(S2);

figure(2);clf;
plot(Q.tVec,Q.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 

figure(3);clf;
psiError = unwrap(n*Q.tVec + pi - Q.state.eMat(:,3));
meanPsiErrorInt = round(mean(psiError)/2/pi);
plot(Q.tVec,psiError - meanPsiErrorInt*2*pi);
grid on;
xlabel('Time (sec)');
ylabel('\Delta \psi (rad)');
title('Yaw angle error');

figure(5);clf;
plot(Q.state.rMat(:,1), Q.state.rMat(:,2)); 
axis equal; grid on; hold on;
xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');
for k=1:100:length(Q.state.eMat)
    RBI = euler2dcm(Q.state.eMat(k,:)');
    xI = RBI'*[1;0;0];
    vectorX = [Q.state.rMat(k,1) Q.state.rMat(k,1)+xI(1)];
    vectorY = [Q.state.rMat(k,2) Q.state.rMat(k,2)+xI(2)];
    plot(vectorX, vectorY, '-r.');
end
legend('Horizontal Position', 'Body x axis direction')

figure(6);clf;
errorMat = Q.state.rMat(1:2:end,:)-R.rIstar;
error = zeros(size(errorMat,1),1);
for k = 1:size(error,1)
    error(k) = norm(errorMat(k,:));
end
plot(Q.tVec(1:2:end), error);
xlabel('Time (sec)');
ylabel('Magnitude of error vector (m)');
title('Position Error over Time');


