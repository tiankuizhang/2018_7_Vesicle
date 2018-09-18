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
		INSTANCE = fullfile(HOME,[DATE,'_SurfaceDiffusion_',HOUR,num2str(floor(rand()*1000))]);
	end
	% crete folder for this simulation instance
	mkdir(INSTANCE)
		
% now make subfolders to record src,img,videos etc
	JPG = fullfile(INSTANCE,'jpg'); % store animation frame
	PNG = fullfile(INSTANCE,'PNG'); % store animation frame with transparent background 
	GIF = fullfile(INSTANCE,'GIF'); % animation with transparent background
	VIDEO = fullfile(INSTANCE,'videos');
	MAT = fullfile(INSTANCE,'mat');
	SRC = fullfile(INSTANCE,'src');

	mkdir(JPG);
	mkdir(PNG);
	mkdir(GIF);
	mkdir(VIDEO);
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

%map.plotSurface(0,1,'green',1)

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

textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

for i=1:20000
	
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

	if mod(i,20)==0
		timeStr = [sprintf('%05d', i), ': ', num2str(time)];

		clf

		subplot(1,2,1)
		map.plot	
		ax = gca;
		ax.Visible = 'off';
		text(textX, textY, textZ, timeStr, 'Color', 'red', 'FontSize', 14);

		subplot(1,2,2)
		map.plotSurface(0,1,'Green',1)	
		ax = gca;
		ax.Visible = 'off';
		text(textX, textY, textZ, timeStr, 'Color', 'red', 'FontSize', 14);

		drawnow

		saveas(gcf, fullfile(JPG, [sprintf('%05d',i),'isosurface','.jpg']))
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
% create video from imges

imageNames = dir(fullfile(JPG,'*.jpg'));
imageNames = {imageNames.name}';

videoOutput = fullfile(GIF,'surfaceDiffusion.gif');

for ii = 1:length(imageNames)
	img = imread(fullfile(JPG,imageNames{ii}));

	% save image with transparent background
	alphaChannel = all(img>150,3);
	imwrite(img, fullfile(PNG, [sprintf('%05d',ii), '.png']), 'Alpha', double(~alphaChannel));
	
	% create gif with transparent background
	[A,map] = rgb2ind(img,256);
	BGColor = double(A(1)); % background color to be set to be transparent
	if ii == 1
		imwrite(A, map, videoOutput, 'gif', 'LoopCount', Inf, 'DelayTime', 1, 'TransparentColor', BGColor);
	else
		imwrite(A, map, videoOutput, 'gif', 'WriteMode', 'append', 'DelayTime', 1, 'TransparentColor', BGColor);
	end

end

