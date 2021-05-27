% draw convergence pictures

load('CFL1.mat')
load('CFL0_5.mat')
load('CFL0_25.mat')
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
plot(t0_25(t0_25<time),ene0_25(t0_25<time),'g')
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
plot(t0_25(t0_25<time),P0_25(t0_25<time),'g')
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
plot(t0_25(t0_25<time),ten0_25(t0_25<time),'g')
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
semilogy(t0_25(t0_25<time),abs(da0_25(t0_25<time)),'g-.')
semilogy(t0_1(t0_1<time),abs(da0_1(t0_1<time)),'b-.')
%semilogy(t0_01(t0_01<time),abs(da0_01(t0_01<time)),'g-.')

semilogy(t1(t1<time),abs(dv1(t1<time)),'m--')
semilogy(t0_5(t0_5<time),abs(dv0_5(t0_5<time)),'r--')
semilogy(t0_25(t0_25<time),abs(dv0_25(t0_25<time)),'g--')
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
plot(t0_25(t0_25<time & t0_25>time2),ene0_25(t0_25<time & t0_25>time2),'g')
plot(t0_1(t0_1<time & t0_1>time2),ene0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),ene0_01(t0_01<time & t0_01>time2),'g')
axis tight

axes('position', [.7 .7 .12 .12])
box on
hold on
plot(t1(t1<time & t1>time2),P1(t1<time & t1>time2),'m')
plot(t0_5(t0_5<time & t0_5>time2),P0_5(t0_5<time & t0_5>time2),'r')
plot(t0_25(t0_25<time & t0_25>time2),P0_25(t0_25<time & t0_25>time2),'g')
plot(t0_1(t0_1<time & t0_1>time2),P0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),P0_01(t0_01<time & t0_01>time2),'g.')
axis tight

axes('position', [.3 .2 .12 .12])
box on
hold on
plot(t1(t1<time & t1>time2),ten1(t1<time & t1>time2),'m')
plot(t0_5(t0_5<time & t0_5>time2),ten0_5(t0_5<time & t0_5>time2),'r')
plot(t0_25(t0_25<time & t0_25>time2),ten0_25(t0_25<time & t0_25>time2),'g')
plot(t0_1(t0_1<time & t0_1>time2),ten0_1(t0_1<time & t0_1>time2),'b')
%plot(t0_01(t0_01<time & t0_01>time2),ten0_01(t0_01<time & t0_01>time2),'g')
axis tight

