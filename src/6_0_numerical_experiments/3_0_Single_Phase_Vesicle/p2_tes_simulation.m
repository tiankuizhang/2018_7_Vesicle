% test if NE works

Beta = 0.8;
CFLNumber = 1;
iteration = 2000;
SampleRate = 10;

simu = SD.Simulation(mfilename,	...
		['SurfaceDiffusion:semi-implicit', ...
		 '_Beta_', num2str(Beta), ...
		 '_CFLNumber_', num2str(CFLNumber), ...
		 '_iteration_', num2str(iteration), ...
		 '_SampleRate_', num2str(SampleRate)]);

simu.simulationStart

pwd
ls

SD.NE.SurfaceDiffusion.PerturbedCylinderPinch(Beta,CFLNumber,iteration,SampleRate,true,simu)

simu.simulationEnd

