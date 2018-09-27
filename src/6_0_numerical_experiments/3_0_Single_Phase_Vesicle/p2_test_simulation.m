% test if NE works

Beta = 0.9;
CFLNumber = 1;
iteration = 200;
SampleRate = 1;
ReinitializationRate = 10;

simu = SD.Simulation(mfilename,	...
		['SurfaceDiffusion', ...
		 '_Beta_', num2str(Beta), ...
		 '_CFLNumber_', num2str(CFLNumber), ...
		 '_iteration_', num2str(iteration), ...
		 '_SampleRate_', num2str(SampleRate), ...
		 '_reiniRate_', num2str(ReinitializationRate)]);

simu.simulationStart

pwd
ls

SD.NE.SurfaceDiffusion.PerturbedCylinderPinch(Beta,CFLNumber,iteration,SampleRate,...
		ReinitializationRate,true,simu)

simu.simulationEnd
SD.NE.processImage(30, 'SurfaceDiffusion')

