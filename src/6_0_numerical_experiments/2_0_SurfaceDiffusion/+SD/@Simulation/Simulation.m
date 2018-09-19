%Simulation : prepare a unique environment for current simulation by
% copying all necessary files to a new location
classdef Simulation < handle

	properties
		FILE % name of the m file 
		simulationName % a string to create a unique folder

		INSTANCE % folder to save all simulation
		SRC % folder to save a copy of the source code

		MAT
		JPG
		PNG
		GIF
	end


	methods
		function obj = Simulation(FILE, simulationName)
			obj.FILE = [FILE, '.m'];
			obj.simulationName = simulationName;
		end

		simulationStart(obj)
		simulationEnd(obj)
		processImage(obj, videoTime)

	end

end
