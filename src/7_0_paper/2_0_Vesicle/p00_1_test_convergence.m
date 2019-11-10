% draw convergence pictures

load('oblate64.mat')
load('oblate96.mat')
load('oblate128.mat')

fprintf('grid size 64^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure64, Tension64)
fprintf('grid size 96^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure96, Tension96)
fprintf('grid size 128^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure128, Tension128)

figure
clf
hold on

time = 1e-3;

plot(t64(t64<time),ene64(t64<time))
plot(t96(t96<time),ene96(t96<time))
plot(t128(t128<time),ene128(t128<time))

box on
legend('Grid Size: 64X64X64', 'Grid Size: 96X96X96', 'Grid Size: 128X128X128')
xlabel('time')
ylabel('elastic bending energy')
title('energy vs time for different grid sizes')

