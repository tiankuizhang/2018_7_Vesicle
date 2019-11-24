% draw convergence pictures

load('CFL1.mat')
load('CFL0_5.mat')
load('CFL0_1.mat')
%load('CFL1.mat')
%load('CFL0_01.mat')

figure
clf

time = 0.7e-3;
time2 = 0.5e-3;

subplot(2,2,1)
hold on
plot(t1(t1<time),ene1(t1<time),'m')
plot(t0_5(t0_5<time),ene0_5(t0_5<time),'r')
plot(t0_1(t0_1<time),ene0_1(t0_1<time),'b')
%plot(t0_01(t0_01<time),ene0_01(t0_01<time),'g')

box on
xlabel('time')
ylabel('elastic bending energy')
title('energy vs time for different CFL numbers')


subplot(2,2,2)
hold on
plot(t1(t1<time),P1(t1<time),'m')
plot(t0_5(t0_5<time),P0_5(t0_5<time),'r')
plot(t0_1(t0_1<time),P0_1(t0_1<time),'b')
%plot(t0_01(t0_01<time),P0_01(t0_01<time),'g')

box on
xlabel('time')
ylabel('P')
title('P vs time for different CFL numbers')

subplot(2,2,3)
hold on
plot(t1(t1<time),ten1(t1<time),'m')
plot(t0_5(t0_5<time),ten0_5(t0_5<time),'r')
plot(t0_1(t0_1<time),ten0_1(t0_1<time),'b')
%plot(t0_01(t0_01<time),ten0_01(t0_01<time),'g')

box on
xlabel('time')
ylabel('\sigma')
title('\sigma vs time for different CFL numbers')

subplot(2,2,4)
semilogy(t1(t1<time),abs(da1(t1<time)),'m-.')
hold on
semilogy(t0_5(t0_5<time),abs(da0_5(t0_5<time)),'r-.')
semilogy(t0_1(t0_1<time),abs(da0_1(t0_1<time)),'b-.')
%semilogy(t0_01(t0_01<time),abs(da0_01(t0_01<time)),'g-.')

semilogy(t1(t1<time),abs(dv1(t1<time)),'m--')
semilogy(t0_5(t0_5<time),abs(dv0_5(t0_5<time)),'r--')
semilogy(t0_1(t0_1<time),abs(dv0_1(t0_1<time)),'b--')
%semilogy(t0_01(t0_01<time),abs(dv0_01(t0_01<time)),'g--')

box on
xlabel('time')
ylabel('relative error')
title('e_A and e_V vs time for different CFL numbers')

axes('position', [.3 .7 .12 .12])
box on
hold on
plot(t1(t1<time & t1>time2),ene1(t1<time & t1>time2),'m')
plot(t0_5(t0_5<time & t0_5>time2),ene0_5(t0_5<time & t0_5>time2),'r')
plot(t0_1(t0_1<time & t0_1>time2),ene0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),ene0_01(t0_01<time & t0_01>time2),'g')
axis tight

axes('position', [.7 .7 .12 .12])
box on
hold on
plot(t1(t1<time & t1>time2),P1(t1<time & t1>time2),'m')
plot(t0_5(t0_5<time & t0_5>time2),P0_5(t0_5<time & t0_5>time2),'r')
plot(t0_1(t0_1<time & t0_1>time2),P0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),P0_01(t0_01<time & t0_01>time2),'g.')
axis tight

axes('position', [.3 .2 .12 .12])
box on
hold on
plot(t1(t1<time & t1>time2),ten1(t1<time & t1>time2),'m')
plot(t0_5(t0_5<time & t0_5>time2),ten0_5(t0_5<time & t0_5>time2),'r')
plot(t0_1(t0_1<time & t0_1>time2),ten0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),ten0_01(t0_01<time & t0_01>time2),'g')
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw convergence pictures
%
%load('oblate48.mat')
%load('oblate64.mat')
%load('oblate80.mat')
%load('oblate96.mat')
%
%figure
%clf
%
%time = 1.6e-3;
%time2 = 0.6e-3;
%
%subplot(2,2,1)
%hold on
%plot(t48(t48<time),ene48(t48<time),'m')
%plot(t64(t64<time),ene64(t64<time),'r')
%plot(t80(t80<time),ene80(t80<time),'b')
%plot(t96(t96<time),ene96(t96<time),'g')
%
%box on
%xlabel('time')
%ylabel('elastic bending energy')
%title('energy vs time for different grid sizes')
%
%
%subplot(2,2,2)
%hold on
%plot(t48(t48<time),P48(t48<time),'m')
%plot(t64(t64<time),P64(t64<time),'r')
%plot(t80(t80<time),P80(t80<time),'b')
%plot(t96(t96<time),P96(t96<time),'g')
%
%box on
%xlabel('time')
%ylabel('P')
%title('P vs time for different grid sizes')
%
%subplot(2,2,3)
%hold on
%plot(t48(t48<time),ten48(t48<time),'m')
%plot(t64(t64<time),ten64(t64<time),'r')
%plot(t80(t80<time),ten80(t80<time),'b')
%plot(t96(t96<time),ten96(t96<time),'g')
%
%box on
%xlabel('time')
%ylabel('\sigma')
%title('\sigma vs time for different grid sizes')
%
%subplot(2,2,4)
%semilogy(t48(t48<time),abs(da48(t48<time)),'m-.')
%hold on
%semilogy(t64(t64<time),abs(da64(t64<time)),'r-.')
%semilogy(t80(t80<time),abs(da80(t80<time)),'b-.')
%semilogy(t96(t96<time),abs(da96(t96<time)),'g-.')
%
%semilogy(t48(t48<time),abs(dv48(t48<time)),'m--')
%semilogy(t64(t64<time),abs(dv64(t64<time)),'r--')
%semilogy(t80(t80<time),abs(dv80(t80<time)),'b--')
%semilogy(t96(t96<time),abs(dv96(t96<time)),'g--')
%
%box on
%xlabel('time')
%ylabel('relative error')
%title('e_A and e_V vs time for different grid sizes')
%
%axes('position', [.3 .7 .12 .12])
%box on
%hold on
%plot(t48(t48<time & t48>time2),ene48(t48<time & t48>time2),'m')
%plot(t64(t64<time & t64>time2),ene64(t64<time & t64>time2),'r')
%plot(t80(t80<time & t80>time2),ene80(t80<time & t80>time2),'b')
%plot(t96(t96<time & t96>time2),ene96(t96<time & t96>time2),'g')
%axis tight
%
%axes('position', [.7 .7 .12 .12])
%box on
%hold on
%plot(t48(t48<time & t48>time2),P48(t48<time & t48>time2),'m')
%plot(t64(t64<time & t64>time2),P64(t64<time & t64>time2),'r')
%plot(t80(t80<time & t80>time2),P80(t80<time & t80>time2),'b')
%plot(t96(t96<time & t96>time2),P96(t96<time & t96>time2),'g.')
%axis tight
%
%axes('position', [.3 .2 .12 .12])
%box on
%hold on
%plot(t48(t48<time & t48>time2),ten48(t48<time & t48>time2),'m')
%plot(t64(t64<time & t64>time2),ten64(t64<time & t64>time2),'r')
%plot(t80(t80<time & t80>time2),ten80(t80<time & t80>time2),'b')
%plot(t96(t96<time & t96>time2),ten96(t96<time & t96>time2),'g')
%axis tight
