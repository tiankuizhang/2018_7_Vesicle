clear; clc;

load('CFL0_1')
load('CFL0_25')
load('CFL0_5')
load('CFL1')


startTime = 0.001;
endTime = 0.0016;

mask1 = (t1 > startTime) & (t1 < endTime);
ene1 = ene1(mask1);
P1 = P1(mask1);
ten1 = ten1(mask1);
t1 = t1(mask1);

mask0_5 = (t0_5 > startTime) & (t0_5 < endTime);
ene0_5 = ene0_5(mask0_5);
P0_5 = P0_5(mask0_5);
ten0_5 = ten0_5(mask0_5);
t0_5 = t0_5(mask0_5);

mask0_25 = (t0_25 > startTime) & (t0_25 < endTime);
ene0_25 = ene0_25(mask0_25);
P0_25 = P0_25(mask0_25);
ten0_25 = ten0_25(mask0_25);
t0_25 = t0_25(mask0_25);

mask0_1 = (t0_1 > startTime) & (t0_1 < endTime);
ene0_1 = ene0_1(mask0_1);
P0_1 = P0_1(mask0_1);
ten0_1 = ten0_1(mask0_1);
t0_1 = t0_1(mask0_1);


totalTime1 = t1(end) - t1(1);
t47 = t1(2:end) - t1(1:end-1);
ene47 = (ene1(2:end) + ene1(1:end-1)) * 0.5;
P47 = (P1(2:end) + P1(1:end-1)) * 0.5;
ten47 = (ten1(2:end) + ten1(1:end-1)) * 0.5;
energyTimeIntegral47 = t47*ene47 / totalTime1;
pressureTimeIntegral47 = t47*P47' / totalTime1;
tensionTimeIntegral47 = t47*ten47' / totalTime1;

totalTime0_5 = t0_5(end) - t0_5(1);
t63 = t0_5(2:end) - t0_5(1:end-1);
ene63 = (ene0_5(2:end) + ene0_5(1:end-1)) * 0.5;
P63 = (P0_5(2:end) + P0_5(1:end-1)) * 0.5;
ten63 = (ten0_5(2:end) + ten0_5(1:end-1)) * 0.5;
energyTimeIntegral63 = t63*ene63 / totalTime0_5;
pressureTimeIntegral63 = t63*P63' / totalTime0_5;
tensionTimeIntegral63 = t63*ten63' / totalTime0_5;

totalTime0_25 = t0_25(end) - t0_25(1);
t95 = t0_25(2:end) - t0_25(1:end-1);
ene95 = (ene0_25(2:end) + ene0_25(1:end-1)) * 0.5;
P95 = (P0_25(2:end) + P0_25(1:end-1)) * 0.5;
ten95 = (ten0_25(2:end) + ten0_25(1:end-1)) * 0.5;
energyTimeIntegral95 = t95*ene95 / totalTime0_25;
pressureTimeIntegral95 = t95*P95' / totalTime0_25;
tensionTimeIntegral95 = t95*ten95' / totalTime0_25;


totalTime0_1 = t0_1(end) - t0_1(1);
t27 = t0_1(2:end) - t0_1(1:end-1);
ene27 = (ene0_1(2:end) + ene0_1(1:end-1)) * 0.5;
P27 = (P0_1(2:end) + P0_1(1:end-1)) * 0.5;
ten27 = (ten0_1(2:end) + ten0_1(1:end-1)) * 0.5;
energyTimeIntegral27 = t27*ene27 / totalTime0_1;
pressureTimeIntegral27 = t27*P27' / totalTime0_1;
tensionTimeIntegral27 = t27*ten27' / totalTime0_1;

EER1 = abs(energyTimeIntegral47 - energyTimeIntegral27);
EER0_5 = abs(energyTimeIntegral63 - energyTimeIntegral27);
EER0_25 = abs(energyTimeIntegral95 - energyTimeIntegral27);

PER1 = abs(pressureTimeIntegral47 - pressureTimeIntegral27);
PER0_5 = abs(pressureTimeIntegral63 - pressureTimeIntegral27);
PER0_25 = abs(pressureTimeIntegral95 - pressureTimeIntegral27);

TER1 = abs(tensionTimeIntegral47 - tensionTimeIntegral27);
TER0_5 = abs(tensionTimeIntegral63 - tensionTimeIntegral27);
TER0_25 = abs(tensionTimeIntegral95 - tensionTimeIntegral27);

fprintf('\n')
fprintf('\t CFL number \tenergy equilibrium \tpressure equilibrium \ttension equilibrium\n');
fprintf('\t %5.2f:   \t\t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 1, energyTimeIntegral47, pressureTimeIntegral47, tensionTimeIntegral47);
fprintf('\t %5.2f:   \t\t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 0.5, energyTimeIntegral63, pressureTimeIntegral63, tensionTimeIntegral63);
fprintf('\t %5.2f:   \t\t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 0.25, energyTimeIntegral95, pressureTimeIntegral95, tensionTimeIntegral95);
fprintf('\t %5.2f:   \t\t%10.5e,  \t\t%10.5e, \t\t%10.5e\n',0.1, energyTimeIntegral27, pressureTimeIntegral27, tensionTimeIntegral27);

fprintf('\n')
fprintf('\t CFL number\tenergy error \t order \t\tpressure error \t order\ttension error \t order\n');
fprintf('\t %5.2f \t\t %5.3e \t -     \t\t %5.3e \t - \t %5.3e \t - \n',1,EER1,PER1,TER1);
fprintf('\t %5.2f \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 0.5,EER0_5,log(EER1/EER0_5)/log(2),PER0_5,log(PER1/PER0_5)/log(2), TER0_5,log(TER1/TER0_5)/log(2));
fprintf('\t %5.2f \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 0.25,EER0_25,log(EER0_5/EER0_25)/log(2),PER0_25,log(PER0_5/PER0_25)/log(2), TER0_25,log(TER0_5/TER0_25)/log(2));

fid = fopen('JCP_temperal_convergence_test.txt',  'w')

fprintf(fid, '\t CFL number \tenergy equilibrium \tpressure equilibrium \ttension equilibrium\n');
fprintf(fid, '\t %5.2f:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 1, energyTimeIntegral47, pressureTimeIntegral47, tensionTimeIntegral47);
fprintf(fid, '\t %5.2f:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 0.5, energyTimeIntegral63, pressureTimeIntegral63, tensionTimeIntegral63);
fprintf(fid, '\t %5.2f:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 0.25, energyTimeIntegral95, pressureTimeIntegral95, tensionTimeIntegral95);
fprintf(fid, '\t %5.2f:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n',0.1, energyTimeIntegral27, pressureTimeIntegral27, tensionTimeIntegral27);

fprintf(fid, '\n');
fprintf(fid, '\t CFL number \tenergy error \t order \t\tpressure error \t order\ttension error \t order\n');
fprintf(fid, '\t %5.2f \t\t %5.3e \t -     \t\t %5.3e \t - \t %5.3e \t - \n',1,EER1,PER1,TER1);
fprintf(fid, '\t %5.2f \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 0.5,EER0_5,log(EER1/EER0_5)/log(2),PER0_5,log(PER1/PER0_5)/log(2), TER0_5,log(TER1/TER0_5)/log(2));
fprintf(fid, '\t %5.2f \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 0.25,EER0_25,log(EER0_5/EER0_25)/log(2),PER0_25,log(PER0_5/PER0_25)/log(2), TER0_25,log(TER0_5/TER0_25)/log(2));

fclose(fid);
