% simulate single phase vesicle without imcompressibility

% create the initial distance map
[x,y,z,f] = SD.Shape.Ellipsoid([128,128,64],0.7,"o");
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,f);
map.F = map.WENORK3Reinitialization(map.F,100);

map.GPUsetCalculusToolBox
InitialArea = map.calArea;
InitialVolume = map.calVolume;
ReduceVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);

fprintf('initial area: %4.5f, initial volume: %4.5f, reduced volume: %4.5f\n', ...
		InitialArea, InitialVolume, ReduceVolume)

% name and size of figure
%FIG = figure('Name','Surface Diffusion of perturbed cylinder','Position',[10 10 1600 800])

% position to show iteration number and accumulated time
textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

MaxResolvedCurvature = 2.0 / map.GD3.Ds;

% bending mudulus
Kappa = 1.0;
CFLNumber = 1.0;

% dynamics
time = 0;
for ii = 1:1
	map.GPUsetCalculusToolBox
	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = (CurrentVolume - InitialVolume) / InitialVolume;
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);

	fprintf('area percent error: %4.5f, volume percent error: %4.5f, reduced volume: %4.5f\n', ...
		DiffArea, DiffVolume, ReduceVolume)

	% extend mean curvature away from surface
	mask = abs(map.F)<2*map.GD3.Ds;
	map.MeanCurvature = sign(map.MeanCurvature) .* ...
		min(MaxResolvedCurvature, abs(map.MeanCurvature));
	MaxCurvatureBeforeExtend = max(abs(map.MeanCurvature(mask)));
	
	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	MaxCurvature = max(abs(MeanCurvature(mask)));
	
	% surface Laplacian of mean curvature and numerical Hamiltonian
	MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian(MeanCurvature); 
	NormalSpeedBend = Kappa * (MeanCurvatureSurfaceLaplacian + 0.5 * MeanCurvature .* ...
			(MeanCurvature.^2 - 4 * map.GaussianCurvature) );
	NormalSpeedBend = map.WENORK3Extend(NormalSpeedBend, 100);

	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;


end
