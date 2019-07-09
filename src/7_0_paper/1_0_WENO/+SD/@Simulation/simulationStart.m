function simulationStart(obj)

	% record simulation start time
		TIME = datetime('now')
		DATE = datestr(TIME,'yy_mm_dd');
		HOUR = datestr(TIME, 'HH_MM_SS');

	% make folder to store src,img,video etc.
		% folder to save all simulation instances
		HOME = '/extra/tiankuizhang'; % folder to save all simulation instances
	
		% subfolder to save current simulation, tagged with HOUR and a random number
		obj.INSTANCE = fullfile(HOME, ...
				[DATE, '_', obj.simulationName, '_', HOUR, num2str(floor(rand()*1000))]);
		while exist(obj.INSTANCE) % if INSTANCE folder aready exists, rename it
			obj.INSTANCE = fullfile(HOME, ...
					[DATE, '_', obj.simulationName, '_', HOUR, num2str(floor(rand()*1000))]);
		end
		% crete folder for this simulation instance
		mkdir(obj.INSTANCE)
		
	% make a subfolder to backup all source files
		obj.SRC = fullfile(obj.INSTANCE,'src');

		CLASSDestination = fullfile(obj.SRC, '+SD');
		mkdir(CLASSDestination)
		CLASSFILE = fullfile(pwd, '+SD');
		copyfile(CLASSFILE, CLASSDestination);

		CUDADestination = fullfile(obj.SRC, 'CUDA_Code');
		mkdir(CUDADestination)
		CUDAFILE = fullfile(pwd, 'CUDA_Code'); 
		copyfile(CUDAFILE, CUDADestination);

		CurrentFILE = fullfile(pwd, obj.FILE);
		copyfile(CurrentFILE, obj.SRC);

		cd(obj.SRC)

	% make folders to store itermediate simulation results
		obj.MAT = fullfile(obj.INSTANCE,'mat');
		obj.JPG = fullfile(obj.INSTANCE,'jpg'); 
		obj.PNG = fullfile(obj.INSTANCE,'png');  
		obj.GIF = fullfile(obj.INSTANCE,'gif'); 

		mkdir(obj.JPG);
		mkdir(obj.PNG);
		mkdir(obj.GIF);
		mkdir(obj.MAT);

	% record test instance information
		TestInfo = fopen(fullfile(obj.INSTANCE,'TestStartInfo'), 'w');
		fprintf(TestInfo, 'test start time: ');
		fprintf(TestInfo, [datestr(TIME, 'yy/mm/dd HH:MM:SS'),'\n']);
		fprintf(TestInfo, ['test file: ', obj.FILE, '\n']);
		fclose(TestInfo);
	
	% diary the command window
		diary(fullfile(obj.INSTANCE, 'command_window'))
		diary on
end










