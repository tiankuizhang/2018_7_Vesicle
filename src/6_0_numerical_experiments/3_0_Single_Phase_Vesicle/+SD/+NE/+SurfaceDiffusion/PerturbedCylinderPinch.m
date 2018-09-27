% numerical experiment handy for tests
% this will run a simulation for pinch of a perturbed cylinder
% Archive is a bool class determining wether to archive this simulation instance or not
% if Archive is true, then a simu object is used to save figures etc
function PerturbedCylinderPinch(Beta, CFLNumber, iteration, SampleRate, Archive, simu)

	% create the initial distance map
	[x,y,z,f] = SD.Shape.PerturbedCylinder([128,128,128], Beta);
	grid = SD.GD3(x,y,z);
	map = SD.SDF3(grid,x,y,z,f);
	map.F = map.WENORK3Reinitialization(map.F,100);

	% name and size of figure
	FIG = figure('Name','Surface Diffusion of perturbed cylinder','Position',[10 10 1600 800])

	% position to show iteration number and accumulated time
	textX = gather(map.GD3.xmin);
	textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
	textZ = gather(map.GD3.zmin);

	MaxResolvedCurvature = 2.0 / map.GD3.Ds;

	% dynamics
	time = 0;
	for i=1:iteration
		
		map.GPUsetCalculusToolBox

		% extend mean curvature away from surface
		mask = abs(map.F)<2*map.GD3.Ds;
		map.MeanCurvature = sign(map.MeanCurvature) .* ...
			min(MaxResolvedCurvature, abs(map.MeanCurvature));
		MaxCurvatureBeforeExtend = max(abs(map.MeanCurvature(mask)));
	
		MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
		MaxCurvature = max(abs(MeanCurvature(mask)));
	
		% surface Laplacian of mean curvature and numerical Hamiltonian
		MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian(MeanCurvature); 
		NormalSpeed = - map.WENORK3Extend(MeanCurvatureSurfaceLaplacian .* map.FGradMag,100);
		
		MaxSpeed = max(abs(NormalSpeed(mask)));
		
		Dt = CFLNumber * map.GD3.Ds / MaxSpeed;
		time = time + Dt;
		
		fprintf('iter: %d, time: %4.10f, MaxCrv: %4.5f, MaxCrv: %4.5f Dt: %4.5f\n', i,time, ...
				MaxCurvatureBeforeExtend/MaxResolvedCurvature, ...
				MaxCurvature/MaxResolvedCurvature,Dt/map.GD3.Ds^4)
	
		% now smooth NormalSpeed
		%normalSpeedSmoothed = smoothFFT(map, NormalSpeed, Dt, 0.5 );
		normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, 0.5);

		normalSpeedSmoothed = map.WENORK3Extend(normalSpeedSmoothed, 100);
		map.F = map.F - Dt * normalSpeedSmoothed;
	
		if mod(i,SampleRate)==0 
			timeStr = [sprintf('%04d: %0.5e', i,time)];
	
			clf(FIG)
	
			subplot(1,2,1)
			map.plot	
			ax = gca;
			ax.Visible = 'off';
			th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
			set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')
	
			subplot(1,2,2)
			map.plotSurface(0,1,'Green','none')	
			ax = gca;
			ax.Visible = 'off';
			th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
			set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')
	
			drawnow
	
			if Archive
				saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
			end
		end

	end


end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with FFT
function normalSpeedSmoothed = smoothFFT(map, NormalSpeed, Dt, Alpha)

	normalSpeedFFT = fftn(NormalSpeed);
	normalSpeedFFTsmoothed = normalSpeedFFT ./ ( 1 + Alpha * Dt *	...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
	normalSpeedSmoothed = real(ifftn(normalSpeedFFTsmoothed));

end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt + Alpha * Dt * map.GD3.LBiLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

%	normalSpeedSmoothed = gmres(@afun, S, [], 1e-12, 300, @mfun);
	normalSpeedSmoothed = gmres(Op, S, [], 1e-12, 300, @mfun);
	normalSpeedSmoothed = reshape(normalSpeedSmoothed, map.GD3.Size);

	% functional form for Op
	function y = afun(x)
		y = Op * x;
	end
	% preconditioner
	function y = mfun(S)
		fftS = fftn(reshape(S,map.GD3.Size));
		fftS = fftS ./ (1 + Alpha * Dt * ...
				(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
		y = real(ifftn(fftS));
		y = reshape(y, [map.GD3.NumElt, 1]);
	end

end
