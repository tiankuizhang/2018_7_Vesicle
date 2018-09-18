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

map.simulationStart('testSimulation');

map.F = map.WENORK3Reinitialization(map.F,100);

% name and size of figure
figure('Name','Surface Diffusion of perturbed cylinder','Position',[10 10 1600 800])

time = 0;

MaxResolvedCurvature = 2.0 / map.GD3.Ds;

textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

for i=1:5
	
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
		timeStr = [sprintf('%04d: %0.5e', i,time)];

		clf

		subplot(1,2,1)
		map.plot	
		ax = gca;
		ax.Visible = 'off';
		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')

		subplot(1,2,2)
		map.plotSurface(0,1,'Green',1)	
		ax = gca;
		ax.Visible = 'off';
		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')

		drawnow

		saveas(gcf, fullfile(map.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
	end

end

map.simulationEnd(10);
