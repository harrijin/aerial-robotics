% Top-level script for calling simulateQuadrotorDynamics
clear; clc;
period = 3; % in seconds
radius = 1;
% Total simulation time, in seconds
Tsim = period;
% Update interval, in seconds.  This value should be small relative to the
% shortest time constant of your system.
delt = 0.005;
% Time vector, in seconds 
N = floor(Tsim/delt);
S.tVec = [0:N-1]'*delt;
% Matrix of disturbance forces acting on the body, in Newtons, expressed in I
S.distMat = zeros(N-1,3);
% Initial position in m
S.state0.r = [0 1 0.5]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e = [atan(4*pi^2*radius/period^2/9.8) 0 0]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v = [(2*pi*radius/period) 0 0]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
eDot = [0 0 -1*2*pi*radius/period]';
phi = S.state0.e(1);
theta = S.state0.e(2);
psi = S.state0.e(3);
S.state0.omegaB = [cos(theta) 0 -sin(theta)*cos(phi);...
                   0 1 sin(phi);...
                   sin(theta) 0 cos(theta)*cos(phi)]*eDot;
% Oversampling factor
S.oversampFact = 10;
% Quadrotor parameters and constants
quadParamsScript;
constantsScript;
S.quadParams = quadParams;
S.constants = constants;
% Rotor speeds at each time, in rad/s
kF = S.quadParams.kF;
kN = S.quadParams.kN;
r = S.quadParams.rotor_loc;
J = S.quadParams.Jq;
m = S.quadParams.m;
g = S.constants.g;
omegaB = S.state0.omegaB;
A = [kF(1)*r(2,1) kF(2)*r(2,2) kF(3)*r(2,3) kF(4)*r(2,4);
     -kF(1)*r(1,1) -kF(2)*r(1,2) -kF(3)*r(1,3) -kF(4)*r(1,4);
     -kN' .* S.quadParams.omegaRdir;
     kF(1) kF(2) kF(3) kF(4)];
b = [omegaB(2)*omegaB(3)*(J(3,3)-J(2,2));
     omegaB(1)*omegaB(3)*(J(1,1)-J(3,3));
     omegaB(2)*omegaB(1)*(J(2,2)-J(1,1));
     m*sqrt((2*pi*radius/period)^4+g^2)];
S.omegaMat = (((A\b).^(1/2))*ones(1, N-1))';
% load('Stest.mat');
P = simulateQuadrotorDynamics(S);
% load('Ptest.mat');
S2.tVec = P.tVec;
S2.rMat = P.state.rMat;
S2.eMat = P.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = true;
S2.gifFileName = 'testGif.gif';
S2.bounds=1*[-1 1 -1 1 -1 1];
visualizeQuad(S2);

figure(1);clf;
plot(P.tVec,P.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 

figure(2);clf;
plot(P.state.rMat(:,1), P.state.rMat(:,2)); 
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');