% simulate single phase vesicle without imcompressibility
% using high order schemes 
% it seems that directly using a larger stencil will not give the correct dynamics

% create the initial distance map
[x,y,z,f] = SD.Shape.Ellipsoid([128,128,64],0.65,"o");
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,f);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);

map.GPUsetCalculusToolBox
InitialArea = map.calArea;
InitialVolume = map.calVolume;
ReduceVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);

fprintf('initial area: %4.5f, initial volume: %4.5f, reduced volume: %4.5f\n', ...
		InitialArea, InitialVolume, ReduceVolume)

% name and size of figure
FIG = figure('Name','Single Phase Vesicle','Position',[10 10 1600 800])

% position to show iteration number and accumulated time
textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

MaxResolvedCurvature = 2.0 / map.GD3.Ds;

% bending mudulus
Kappa = 1.0;
CFLNumber = .1;

% dynamics
time = 0;
for i = 1:20000
	map.GPUsetCalculusToolBox
	%map.setCalculusToolBox4
	%map.setCalculusToolBoxWENO
	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - InitialVolume) / InitialVolume;
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);


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

	% now solve for tension and pressure to enfore total area and volume
	c11 = map.surfaceIntegral(MeanCurvature.^2); % surface integral of mean curvature squared
	c12 = map.surfaceIntegral(MeanCurvature); % surface integral of mean curvature
	c21 = c12;
	c22 = InitialArea;

	s1 = - map.surfaceIntegral(MeanCurvature .* NormalSpeedBend) ...
		 - (InitialArea - CurrentArea) / Dt;
	s2 = map.surfaceIntegral(NormalSpeedBend) + (InitialVolume - CurrentVolume) / Dt;
	
	TP = [c11,c12;c21,c22]\[s1;s2];
	Tension = TP(1);
	Pressure = TP(2);

	fprintf('iter: %d, ene: %4.5f, ar: %4.5f, vol: %4.5f, rd: %4.5f\n', ...
			i, c11, DiffArea, DiffVolume, ReduceVolume)
	% now calculate normal Speed
	normalSpeed = (Tension * MeanCurvature - NormalSpeedBend + Pressure) .* map.FGradMag;
	normalSpeed = map.WENORK3Extend(normalSpeed, 100);
	normalSpeedSmoothed = smoothGMRES(map, normalSpeed, Dt, 0.5);
	normalSpeedSmoothed = map.WENORK3Extend(normalSpeedSmoothed, 100);

	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance

	if mod(i,100)==0 
		timeStr = [sprintf('%04d: %0.5e', i,time)];

		clf(FIG)

		map.plotSurface(0,1,'Green','black')	
		ax = gca;
		ax.Visible = 'off';
		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')

		zoom(1.3)

		drawnow

		if false 
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,10)==0
		map.F = map.WENO5RK3Reinitialization(map.F,100);
	end


end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt + Alpha * Dt * map.GD3.LBiLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

	normalSpeedSmoothed = gmres(Op, S, [], 1e-12, 300, @mfun);
	normalSpeedSmoothed = reshape(normalSpeedSmoothed, map.GD3.Size);

	% preconditioner
	function y = mfun(S)
		fftS = fftn(reshape(S,map.GD3.Size));
		fftS = fftS ./ (1 + Alpha * Dt * ...
				(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
		y = real(ifftn(fftS));
		y = reshape(y, [map.GD3.NumElt, 1]);
	end

end
