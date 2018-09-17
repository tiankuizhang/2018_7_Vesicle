%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% run simulation, save imges and create videos

% record simulation start time
	TIME = datetime('now')
	DATE = datestr(TIME,'yy_mm_dd');
	HOUR = datestr(TIME, 'HH_MM_SS');
	FILE = mfilename

% make folder to store src,img,video etc.
	% folder to save all simulation instances
	HOME = '/extra/tiankuizhang'; % folder to save all simulation instances

	% subfolder to save current simulation, tagged with HOUR and a random number
	INSTANCE = fullfile(HOME,[DATE,'_SurfaceDiffusion_',HOUR,num2str(floor(rand()*1000))]);
	while exist(INSTANCE) % if INSTANCE folder aready exists, rename it
		INSTANCE = fullfile(HOME,[DATE,'SurfaceDiffusion',HOUR,num2str(floor(rand()*1000))]);
	end
	% crete folder for this simulation instance
	mkdir(INSTANCE)
		
% now make subfolders to record src,img,videos etc
	IMG = fullfile(INSTANCE,'imges');
	MAT = fullfile(INSTANCE,'mat');
	SRC = fullfile(INSTANCE,'src');

	mkdir(IMG);
	mkdir(MAT);
	mkdir(SRC);

	copyfile(pwd, SRC); % copy source file

% record test instance information
	TestInfo = fopen(fullfile(INSTANCE,'TestStartInfo'), 'w');
	fprintf(TestInfo, 'test start time: ');
	fprintf(TestInfo, [datestr(TIME, 'yy/mm/dd HH:MM:SS'),'\n']);
	fprintf(TestInfo, ['test file: ', FILE, '\n']);
	fclose(TestInfo);

% diary the command window
	diary(fullfile(INSTANCE, 'command_window'))
	diary on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve surface diffusion problem with an explicit method
D = 0.5;
Alpha = 2;

xv = linspace(-2*pi-D,2*pi+D,128)/Alpha;
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

% a cylinder with a sinusoidally varying radius
C1 = 0.5 * pi/Alpha; % maximum radius
C2 = 0.90 * C1/2; % C1 - 2*C2 will be the minimum radius

F1 = sqrt(x.^2+y.^2) - (C1-C2*(cos(Alpha * z) + 1));
F2 = max(z-pi/Alpha,-z-pi/Alpha);

F3 = sqrt(x.^2+y.^2+(z-pi/Alpha).^2) - C1;
F4 = sqrt(x.^2+y.^2+(z+pi/Alpha).^2) - C1;

F5 = min(F3,F4);

F = max(F1,F2);

F = min(F,F5);


% set up the distance function
map = SD.SDF3(grid, x, y, z, F);

map.F = map.WENORK3Reinitialization(map.F,100);

map.plotSurface(0,1,'green',1)

%load('Extend.mat')
%
%grid = SD.GD3(x,y,z);
%map = SD.SDF3(grid,x,y,z,F);
%
%map.GPUsetCalculusToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% explicit scheme

figure(1)

time = 0;

MaxResolvedCurvature = 2.0 / map.GD3.Ds;

for i=1:20
	
	map.GPUsetCalculusToolBox
	mask = abs(map.F)<2*map.GD3.Ds;
	
	map.MeanCurvature = sign(map.MeanCurvature) .* min(MaxResolvedCurvature, abs(map.MeanCurvature));
	MaxCurvatureBeforeExtend = max(abs(map.MeanCurvature(mask)));

	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	%MeanCurvature = map.ENORK2Extend(map.MeanCurvature,100);
	MaxCurvature = max(abs(MeanCurvature(mask)));

	% surface Laplacian of mean curvature
	
	MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian(MeanCurvature); 
	NormalSpeed = - map.WENORK3Extend(MeanCurvatureSurfaceLaplacian .* map.FGradMag,100);
	
	MaxSpeed = max(abs(NormalSpeed(mask)));
	
	CFLNumber = 0.001;
	
	Dt = CFLNumber * map.GD3.Ds / MaxSpeed;
	time = time + Dt;
	
	fprintf('iter: %d, time: %4.10f, MaxCrv: %4.5f, MaxCrv: %4.5f Dt: %4.5f\n', i,time, ...
			MaxCurvatureBeforeExtend/MaxResolvedCurvature, ...
			MaxCurvature/MaxResolvedCurvature,Dt/map.GD3.Ds^4)

	map.F = map.F - Dt * NormalSpeed;
	% there is a sign error in the calculation of MeanCurvature
	% thus here we shall use a plus sign

	if mod(i,1)==0
		clf
		%map.plotSurface(0,1,'Green',1)	
		map.plot	
		title([sprintf('%05d',i), ':', num2str(time)])
		drawnow

		% save figure
		saveas(gcf, fullfile(IMG, [sprintf('%05d',i),'isosurface','.jpg']))
		% save current distance map
		DistanceMap = map.F;
		save(fullfile(MAT,['DistanceMap',sprintf('%05d',i),'.mat']), 'DistanceMap');

	end

	if mod(i,20)==0
		%map.F = map.WENORK3Reinitialization(map.F,100);
	end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log test end information
	TestInfo = fopen(fullfile(INSTANCE,'TestEndInfo'), 'w');
	fprintf(TestInfo, 'test end time: \t');
	fprintf(TestInfo, [datestr(datetime('now'),'yy/mm/dd HH:MM:SS'), '\n']);
	fclose(TestInfo)

	diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

