% prepare for simulation
function simulationStart(obj, simulationName)

	obj.simulationName = simulationName;
		
	% record simulation start time
		TIME = datetime('now')
		DATE = datestr(TIME,'yy_mm_dd');
		HOUR = datestr(TIME, 'HH_MM_SS');
		FILE = mfilename

	% make folder to store src,img,video etc.
		% folder to save all simulation instances
		HOME = '/extra/tiankuizhang'; % folder to save all simulation instances
	
		% subfolder to save current simulation, tagged with HOUR and a random number
		INSTANCE = fullfile(HOME,[DATE, '_', simulationName, '_', HOUR, num2str(floor(rand()*1000))]);
		while exist(INSTANCE) % if INSTANCE folder aready exists, rename it
			INSTANCE = fullfile(HOME,[DATE, '_', simulationName, '_', HOUR, num2str(floor(rand()*1000))]);
		end
		% crete folder for this simulation instance
		mkdir(INSTANCE)
		
		obj.INSTANCE = INSTANCE;
		
	% now make subfolders to record src,img,videos etc
		obj.JPG = fullfile(INSTANCE,'jpg'); % store animation frame
		obj.PNG = fullfile(INSTANCE,'png'); % store animation frame with transparent background 
		obj.GIF = fullfile(INSTANCE,'gif'); % animation with transparent background
		obj.MAT = fullfile(INSTANCE,'mat');
		obj.SRC = fullfile(INSTANCE,'src');
	
		mkdir(obj.JPG);
		mkdir(obj.PNG);
		mkdir(obj.GIF);
		mkdir(obj.MAT);
		mkdir(obj.SRC);
	
		copyfile(pwd, obj.SRC); % copy source file

	% record test instance information
		TestInfo = fopen(fullfile(INSTANCE,'TestStartInfo'), 'w');
		fprintf(TestInfo, 'test start time: ');
		fprintf(TestInfo, [datestr(TIME, 'yy/mm/dd HH:MM:SS'),'\n']);
		fprintf(TestInfo, ['test file: ', FILE, '\n']);
		fclose(TestInfo);
	
	% diary the command window
		diary(fullfile(INSTANCE, 'command_window'))
		diary on
end
