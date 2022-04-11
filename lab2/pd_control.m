kP = 10;
kD = 10;
kY = 1;
P = tf([kY],[1 0 0]);
C = tf([kD kP], [1]);
G = (P*C)/(1+P*C);
step(G)
grid on
stepinfo(G)