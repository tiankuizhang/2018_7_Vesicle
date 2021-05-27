clear; clc;

load('oblate48.mat')
load('oblate64.mat')
load('oblate80.mat')
load('oblate96.mat')
load('oblate128.mat')

startTime = 0.001;
endTime = 0.0016;

mask48 = (t48 > startTime) & (t48 < endTime);
ene48 = ene48(mask48);
P48 = P48(mask48);
ten48 = ten48(mask48);
t48 = t48(mask48);

mask64 = (t64 > startTime) & (t64 < endTime);
ene64 = ene64(mask64);
P64 = P64(mask64);
ten64 = ten64(mask64);
t64 = t64(mask64);

mask80 = (t80 > startTime) & (t80 < endTime);
ene80 = ene80(mask80);
P80 = P80(mask80);
ten80 = ten80(mask80);
t80 = t80(mask80);

mask96 = (t96 > startTime) & (t96 < endTime);
ene96 = ene96(mask96);
P96 = P96(mask96);
ten96 = ten96(mask96);
t96 = t96(mask96);

mask128 = (t128 > startTime) & (t128 < endTime);
ene128 = ene128(mask128);
P128 = P128(mask128);
ten128 = ten128(mask128);
t128 = t128(mask128);


totalTime48 = t48(end) - t48(1);
t47 = t48(2:end) - t48(1:end-1);
ene47 = (ene48(2:end) + ene48(1:end-1)) * 0.5;
P47 = (P48(2:end) + P48(1:end-1)) * 0.5;
ten47 = (ten48(2:end) + ten48(1:end-1)) * 0.5;
energyTimeIntegral47 = t47*ene47 / totalTime48;
pressureTimeIntegral47 = t47*P47' / totalTime48;
tensionTimeIntegral47 = t47*ten47' / totalTime48;

totalTime64 = t64(end) - t64(1);
t63 = t64(2:end) - t64(1:end-1);
ene63 = (ene64(2:end) + ene64(1:end-1)) * 0.5;
P63 = (P64(2:end) + P64(1:end-1)) * 0.5;
ten63 = (ten64(2:end) + ten64(1:end-1)) * 0.5;
energyTimeIntegral63 = t63*ene63 / totalTime64;
pressureTimeIntegral63 = t63*P63' / totalTime64;
tensionTimeIntegral63 = t63*ten63' / totalTime64;

totalTime80 = t80(end) - t80(1);
t79 = t80(2:end) - t80(1:end-1);
ene79 = (ene80(2:end) + ene80(1:end-1)) * 0.5;
P79 = (P80(2:end) + P80(1:end-1)) * 0.5;
ten79 = (ten80(2:end) + ten80(1:end-1)) * 0.5;
energyTimeIntegral79 = t79*ene79 / totalTime80;
pressureTimeIntegral79 = t79*P79' / totalTime80;
tensionTimeIntegral79 = t79*ten79' / totalTime80;

totalTime96 = t96(end) - t96(1);
t95 = t96(2:end) - t96(1:end-1);
ene95 = (ene96(2:end) + ene96(1:end-1)) * 0.5;
P95 = (P96(2:end) + P96(1:end-1)) * 0.5;
ten95 = (ten96(2:end) + ten96(1:end-1)) * 0.5;
energyTimeIntegral95 = t95*ene95 / totalTime96;
pressureTimeIntegral95 = t95*P95' / totalTime96;
tensionTimeIntegral95 = t95*ten95' / totalTime96;


totalTime128 = t128(end) - t128(1);
t27 = t128(2:end) - t128(1:end-1);
ene27 = (ene128(2:end) + ene128(1:end-1)) * 0.5;
P27 = (P128(2:end) + P128(1:end-1)) * 0.5;
ten27 = (ten128(2:end) + ten128(1:end-1)) * 0.5;
energyTimeIntegral27 = t27*ene27 / totalTime128;
pressureTimeIntegral27 = t27*P27' / totalTime128;
tensionTimeIntegral27 = t27*ten27' / totalTime128;

EER48 = abs(energyTimeIntegral47 - energyTimeIntegral27);
EER64 = abs(energyTimeIntegral63 - energyTimeIntegral27);
EER80 = abs(energyTimeIntegral79 - energyTimeIntegral27);
EER96 = abs(energyTimeIntegral95 - energyTimeIntegral27);

PER48 = abs(pressureTimeIntegral47 - pressureTimeIntegral27);
PER64 = abs(pressureTimeIntegral63 - pressureTimeIntegral27);
PER80 = abs(pressureTimeIntegral79 - pressureTimeIntegral27);
PER96 = abs(pressureTimeIntegral95 - pressureTimeIntegral27);

TER48 = abs(tensionTimeIntegral47 - tensionTimeIntegral27);
TER64 = abs(tensionTimeIntegral63 - tensionTimeIntegral27);
TER80 = abs(tensionTimeIntegral79 - tensionTimeIntegral27);
TER96 = abs(tensionTimeIntegral95 - tensionTimeIntegral27);

fprintf('\n')
fprintf('\t gridsize \tenergy equilibrium \tpressure equilibrium \ttension equilibrium\n');
fprintf('\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 48, energyTimeIntegral47, pressureTimeIntegral47, tensionTimeIntegral47);
fprintf('\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 64, energyTimeIntegral63, pressureTimeIntegral63, tensionTimeIntegral63);
fprintf('\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 80, energyTimeIntegral79, pressureTimeIntegral79, tensionTimeIntegral79);
fprintf('\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 96, energyTimeIntegral95, pressureTimeIntegral95, tensionTimeIntegral95);
fprintf('\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n',128, energyTimeIntegral27, pressureTimeIntegral27, tensionTimeIntegral27);

fprintf('\n')
fprintf('\t gridsize \tenergy error \t order \t\tpressure error \t order\ttension error \t order\n');
fprintf('\t %03d \t\t %5.3e \t -     \t\t %5.3e \t - \t %5.3e \t - \n',48,EER48,PER48,TER48)
fprintf('\t %03d \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 64,EER64,log(EER48/EER64)/log(64/48),PER64,log(PER48/PER64)/log(64/48), TER64,log(TER48/TER64)/log(64/48))
fprintf('\t %03d \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 80,EER80,log(EER64/EER80)/log(80/64),PER80,log(PER64/PER80)/log(80/64), TER80,log(TER64/TER80)/log(80/64))
fprintf('\t %03d \t\t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 96,EER96,log(EER80/EER96)/log(96/80),PER96,log(PER80/PER96)/log(96/80), TER96,log(TER80/TER96)/log(96/80))

fid = fopen('JCP_spatial_convergence_test.txt',  'w')

fprintf(fid, '\t gridsize \tenergy equilibrium \tpressure equilibrium \ttension equilibrium\n');
fprintf(fid, '\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 48, energyTimeIntegral47, pressureTimeIntegral47, tensionTimeIntegral47);
fprintf(fid, '\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 64, energyTimeIntegral63, pressureTimeIntegral63, tensionTimeIntegral63);
fprintf(fid, '\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 80, energyTimeIntegral79, pressureTimeIntegral79, tensionTimeIntegral79);
fprintf(fid, '\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n', 96, energyTimeIntegral95, pressureTimeIntegral95, tensionTimeIntegral95);
fprintf(fid, '\t %03d:   \t%10.5e,  \t\t%10.5e, \t\t%10.5e\n',128, energyTimeIntegral27, pressureTimeIntegral27, tensionTimeIntegral27);

fprintf(fid, '\n');
fprintf(fid, '\t gridsize \tenergy error \t order \t\tpressure error \t order\ttension error \t order\n');
fprintf(fid, '\t %03d \t %5.3e \t -     \t\t %5.3e \t - \t %5.3e \t - \n',48,EER48,PER48,TER48);
fprintf(fid, '\t %03d \t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 64,EER64,log(EER48/EER64)/log(64/48),PER64,log(PER48/PER64)/log(64/48), TER64,log(TER48/TER64)/log(64/48));
fprintf(fid, '\t %03d \t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 80,EER80,log(EER64/EER80)/log(80/64),PER80,log(PER64/PER80)/log(80/64), TER80,log(TER64/TER80)/log(80/64));
fprintf(fid, '\t %03d \t %5.3e \t %7.2f \t %5.3e \t %5.2f \t %5.3e \t %5.2f \n', 96,EER96,log(EER80/EER96)/log(96/80),PER96,log(PER80/PER96)/log(96/80), TER96,log(TER80/TER96)/log(96/80));

fclose(fid);
