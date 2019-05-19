function simulationEnd(obj)

	% log test end information
		TestInfo = fopen(fullfile(obj.INSTANCE,'TestEndInfo'), 'w');
		fprintf(TestInfo, 'test end time: \t');
		fprintf(TestInfo, [datestr(datetime('now'),'yy/mm/dd HH:MM:SS'), '\n']);
		fclose(TestInfo)
	
		diary off
end
