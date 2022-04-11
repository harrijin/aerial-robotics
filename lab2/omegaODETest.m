quadParamsScript;
motor = tf([quadParams.cm(1)],[quadParams.taum(1) 1]);
step(motor)
stepinfo(motor)
grid on
%% Ode
omegaDot = @(eA, omega) (quadParams.cm(1)*eA-omega)/quadParams.taum(1);
[t, omegaOut] = ode45(@(t,omega)omegaDot(1,omega), 0:0.001:0.45, 0);
figure(2)
plot(t,omegaOut);
grid on
title('Step Response from ODE');
xlabel('Time (seconds)');
ylabel('Amplitude');