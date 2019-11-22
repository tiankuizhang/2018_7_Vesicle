% draw convergence pictures

load('oblate64.mat')
load('oblate96.mat')
load('oblate128.mat')

%fprintf('grid size 64^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure64, Tension64)
%fprintf('grid size 96^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure96, Tension96)
%fprintf('grid size 128^3, \t pressure: %.2f, \t tension: %.2f \n', Pressure128, Tension128)

figure
clf

time = 3e-3;
time2 = 2.7e-3;

subplot(2,2,1)
hold on
plot(t64(t64<time),ene64(t64<time),'r')
plot(t96(t96<time),ene96(t96<time),'b')
plot(t128(t128<time),ene128(t128<time),'g')

box on
xlabel('time')
ylabel('elastic bending energy')
title('energy vs time for different grid sizes')


subplot(2,2,2)
hold on
plot(t64(t64<time),P64(t64<time),'r')
plot(t96(t96<time),P96(t96<time),'b')
plot(t128(t128<time),P128(t128<time),'g')

box on
xlabel('time')
ylabel('P')
title('P vs time for different grid sizes')

subplot(2,2,3)
hold on
plot(t64(t64<time),ten64(t64<time),'r')
plot(t96(t96<time),ten96(t96<time),'b')
plot(t128(t128<time),ten128(t128<time),'g')

box on
xlabel('time')
ylabel('\sigma')
title('\sigma vs pressure for different grid sizes')

subplot(2,2,4)
%hold on
%plot(t64(t64<time),da64(t64<time),'r')
%plot(t96(t96<time),da96(t96<time),'b')
%plot(t128(t128<time),da128(t128<time),'g')
%plot(t64(t64<time),dv64(t64<time),'r--')
%plot(t96(t96<time),dv96(t96<time),'b--')
%plot(t128(t128<time),dv128(t128<time),'g--')
semilogy(t64(t64<time),abs(da64(t64<time)),'r')
hold on
semilogy(t96(t96<time),abs(da96(t96<time)),'b')
semilogy(t128(t128<time),abs(da128(t128<time)),'g')
semilogy(t64(t64<time),abs(dv64(t64<time)),'r--')
semilogy(t96(t96<time),abs(dv96(t96<time)),'b--')
semilogy(t128(t128<time),abs(dv128(t128<time)),'g--')

box on
xlabel('time')
ylabel('relative error')
title('e_A and e_V vs time for different grid sizes')

axes('position', [.3 .7 .12 .12])
box on
hold on
plot(t64(t64<time & t64>time2),ene64(t64<time & t64>time2),'r')
plot(t96(t96<time & t96>time2),ene96(t96<time & t96>time2),'b')
plot(t128(t128<time & t128>time2),ene128(t128<time & t128>time2),'g')
axis tight

axes('position', [.7 .7 .12 .12])
box on
hold on
plot(t64(t64<time & t64>time2),P64(t64<time & t64>time2),'r')
plot(t96(t96<time & t96>time2),P96(t96<time & t96>time2),'b')
plot(t128(t128<time & t128>time2),P128(t128<time & t128>time2),'g')
axis tight

axes('position', [.3 .2 .12 .12])
box on
hold on
plot(t64(t64<time & t64>time2),ten64(t64<time & t64>time2),'r')
plot(t96(t96<time & t96>time2),ten96(t96<time & t96>time2),'b')
plot(t128(t128<time & t128>time2),ten128(t128<time & t128>time2),'g')
axis tight
