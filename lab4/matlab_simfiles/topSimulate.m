% Top-level script for calling simulateQuadrotorDynamics
clear; clc;
period = 10; % in seconds
radius = 4;
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
S.state0.r = [0 radius 0.5]';
% Initial attitude expressed as Euler angles, in radians
S.state0.e = [0 0 -pi/2]';
% Initial velocity of body with respect to I, expressed in I, in m/s
S.state0.v = [0 0 0]';
% Initial angular rate of body with respect to I, expressed in B, in rad/s
S.state0.omegaB = [0 0 0]';
% Oversampling factor
S.oversampFact = 10;
% Quadrotor parameters and constants
quadParamsScript;
constantsScript;
S.quadParams = quadParams;
S.constants = constants;
% Rotor speeds at each time, in rad/s
% omega0 = 800;
% Non high-fidelity simulator
% S.omegaMat = omega0*ones(4, N-1)';
% P = simulateQuadrotorDynamics(S);
% % High fidelity simulator
% S.eaMat = ((omega0./quadParams.cm)*ones(1, N-1))';
% P = simulateQuadrotorDynamicsHF(S);
% Control simulator
R.tVec = S.tVec;
R.rIstar = [radius*sin(R.tVec*2*pi/period)';radius*cos(R.tVec*2*pi/period)';0.5*ones(1,N)]';
R.vIstar = [2*pi/period*radius*cos(R.tVec*2*pi/period)';2*pi/period*-1*radius*sin(R.tVec*2*pi/period)';zeros(1,N)]';
R.aIstar = [4*pi^2/period^2*-1*radius*sin(R.tVec*2*pi/period)';4*pi^2/period^2*-1*radius*cos(R.tVec*2*pi/period)';zeros(1,N)]';
R.xIstar = -1*[sin(R.tVec*2*pi/period)';cos(R.tVec*2*pi/period)';zeros(1,N)]';
P1.quadParams = quadParams;
P1.constants = constants;
P = simulateQuadrotorControl(R,S,P1);

S2.tVec = P.tVec;
S2.rMat = P.state.rMat;
S2.eMat = P.state.eMat;
S2.plotFrequency = 20;
S2.makeGifFlag = false;
S2.gifFileName = 'testGif.gif';
S2.bounds=1*[-5 5 -5 5 -5 5];
visualizeQuad(S2);

figure(1);clf;
plot(P.tVec,P.state.rMat(:,3)); grid on;
xlabel('Time (sec)');
ylabel('Vertical (m)');
title('Vertical position of CM'); 

figure(2);clf;
axis equal; grid on; hold on;
plot(P.state.rMat(:,1), P.state.rMat(:,2));
% for k=1:400:length(P.state.eMat)
%     RBI = euler2dcm(P.state.eMat(k,:)');
%     xI = RBI'*[1;0;0];
%     vectorX = [P.state.rMat(k,1) P.state.rMat(k,1)+xI(1)];
%     vectorY = [P.state.rMat(k,2) P.state.rMat(k,2)+xI(2)];
%     plot(vectorX, vectorY, '-r.');
% end
% legend('Horizontal Position', 'Body x axis direction')


xlabel('X (m)');
ylabel('Y (m)');
title('Horizontal position of CM');

% figure(3);clf;
% plot(P.tVec,P.state.rMat(:,1)); grid on;
% title('X vs time');
% 
% figure(4);clf;
% plot(P.tVec,P.state.rMat(:,2)); grid on;
% title('Y vs time');
% 
figure(5);clf;
plot(P.tVec,P.state.eMat(:,3)); grid on;
title('Yaw vs time');
xlabel('Time (sec)');
ylabel('Yaw (rad)');

