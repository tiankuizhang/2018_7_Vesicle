% simulate a two phase vesicle conserving total area, volume and area of a particular phase
% incomprssibility is not imposed right now


[x,y,z,f,a,b,c] = SD.Shape.Ellipsoid([96,96,96],0.95,"p");
fprintf('a:%4.5f, b:%4.5f, c:%4.5f\n',a,b,c)
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,f);
map.A = z - 0.70*a;

map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,100,50);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox
InitialArea = map.calArea;
InitialVolume = map.calVolume;
InitialPhaseArea = map.AcalArea;
ReduceVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);

fprintf('initial area: %4.5f, initial volume: %4.5f, reduced volume: %4.5f\n', InitialArea, InitialVolume, ReduceVolume)
% name and size of figure
FIG = figure('Name','MultiPhase Vesicle','Position',[10 10 1600 800])

iso = gather(-20*map.GD3.Ds:4*map.GD3.Ds:20*map.GD3.Ds);

% position to show iteration number and accumulated time
textX = gather(map.GD3.xmin);
%textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textY = gather( map.GD3.ymin );
textZ = gather(map.GD3.zmin);

% bending mudulus
KappaP = 3.0; % bending rigidity for positive psi
KappaN = 1.0; % bending rigidity for negative psi

KappaL = 30; % isotropic line tension
CFLNumber = 0.1;
filterWidth = gather(map.GD3.Ds)*5.0;

% dynamics
time = 0;
for i = 0:3000
	%KappaL = i;
	map.GPUsetCalculusToolBox
	map.GPUAsetCalculusToolBox
	KappaB = KappaP + (KappaN - KappaP) * map.AHeaviside;

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - InitialVolume) / InitialVolume;
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);
	CurrentPhaseArea = map.AcalArea;
	DiffPhaseArea = 100 * (CurrentPhaseArea - InitialPhaseArea) / InitialPhaseArea; 

	% extend mean curvature away from surface
	map.MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);

	% surface Laplacian of mean curvature and numerical Hamiltonian
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(map.MeanCurvature); 
	NormalSpeedBend = KappaB .* (MeanCurvatureSurfaceLaplacian + ...
	 0.5 * map.MeanCurvature .* (map.MeanCurvature.^2 - 4 * map.GaussianCurvature) ) ...
	 - KappaL * map.NormalCurvature .* map.ADiracDelta .* map.AGradMag ;
	NormalSpeedBend = map.ENORK2Extend(NormalSpeedBend, 100);

	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;

	% now solve for tension and pressure to enfore total area and volume
	MeanCurvature = map.MeanCurvature;
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

	% now calculate normal Speed
	normalSpeed = Tension * MeanCurvature - NormalSpeedBend + Pressure;

	% time step level set function
	normalSpeedSmoothed = smoothFFT(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	%normalSpeedSmoothed = smoothGMRES(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);

	%% time step the auxilary level set function
	% determine the Lagrange multiplier that conserves phase area
	GeodesicCurvature = map.AENORK2Extend(map.GeodesicCurvature, 50, 100, 50);
	PhaseLength = map.calLength;
	PhaseChange = (InitialPhaseArea - CurrentPhaseArea) / Dt + ...
		map.AsurfaceIntegral(MeanCurvature .* normalSpeedSmoothed) - ...
		KappaL * map.LineIntegral(map.GeodesicCurvature) ;
	LineTension = - PhaseChange / PhaseLength;

	AnormalSpeed = KappaL * GeodesicCurvature - LineTension;

	AnormalSpeed = smoothFFT(map, AnormalSpeed.*map.AGradMag, Dt, 0.5);
	%AnormalSpeed = smoothGMRES(map, AnormalSpeed.*map.AGradMag, Dt, 0.5);

	% timestep level set function
	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance
	map.A = map.A - Dt * AnormalSpeed;
	map.A = map.WENORK3Extend(map.A, 50);
	map.A = map.WENORK3Reinitialization(map.A, 50);
	map.A = map.WENORK3Extend(map.A, 50);


	c11 = map.surfaceIntegral(KappaB .* MeanCurvature.^2); % surface integral of mean curvature squared
	ene = c11 + KappaL * map.calLength;
	fprintf('iter: %5d, ene: %4.5f, ar: %+4.5f, vol: %+4.5f, rd: %4.5f, pe: %+4.5f\n', ...
			i, ene, DiffArea, DiffVolume, ReduceVolume, DiffPhaseArea)

	if mod(i,20)==0
		timeStr = [sprintf('%04d: %0.5e, %0.5f', i,time,ene)];

		clf(FIG)

		%map.plotSurface(0,1,'Green','black');textZ = gather(map.GD3.zmin);
		%map.plotField(0,normalSpeedSmoothed,0.5)
		map.plotField(0,map.AHeaviside,0.01)
		map.GD3.DrawBox

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')

		zoom(1.0)

		drawnow


		if false 
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,10)==0
		%map.F = map.ENORK2Reinitialization(map.F,100);
		map.F = map.WENO5RK3Reinitialization(map.F,200);
		%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,20,30);
	end

end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt + Alpha * Dt * map.GD3.LBiLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

	[normalSpeedSmoothed,~,~,~,~] = gmres(Op, S, [], 1e-12, 300, @mfun);
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

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothFFT(map, NormalSpeed, Dt, Alpha)

	fftS = fftn(NormalSpeed);
	fftS = fftS ./ (1 + Alpha * Dt * ...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
	normalSpeedSmoothed = real(ifftn(fftS));
	
end




