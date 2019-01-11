% simulate single phase vesicle and measuring local area change 
% impose incompressibility and calculate tension and pressure

% create the initial distance map
%[x,y,z,f] = SD.Shape.Ellipsoid([128,128,64],0.65,"o");
%[x,y,z,f] = SD.Shape.Ellipsoid([64,64,32],0.65,"o");
[x,y,z,f] = SD.Shape.Ellipsoid([64,64,128],0.65,"p");
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
textZ = gather(map.GD3.zmax);

% bending mudulus
Kappa = 1.0;
CFLNumber = 0.1;
filterWidth = gather(map.GD3.Ds)*5.0;
localArea = ones(map.GD3.Size,'gpuArray'); % measurement of local compression and stretch

% dynamics
time = 0;%
for i = 1:100
	map.GPUsetCalculusToolBox

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - InitialVolume) / InitialVolume;
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);

	% extend mean curvature away from surface
	map.MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	
	% surface Laplacian of mean curvature and numerical Hamiltonian
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(map.MeanCurvature); 
	NormalSpeedBend = Kappa * (MeanCurvatureSurfaceLaplacian + 0.5 * map.MeanCurvature ...
			.* (map.MeanCurvature.^2 - 4 * map.GaussianCurvature) );
	NormalSpeedBend = map.ENORK2Extend(NormalSpeedBend, 100);

	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;

	% now solve for local tension and tension and pressure to enfore total area and volume
	c11 = map.surfaceIntegral(map.MeanCurvature.^2); % energy
	c12 = map.surfaceIntegral(map.MeanCurvature); % surface integral of mean curvature
	c21 = c12;
	c22 = InitialArea;
	
	% calculate local tension
	[localTension,residual] = localLagrangeMultiplier(map, -NormalSpeedBend, localArea, Dt, c11/c22); 
	%localTension = localLagrangeMultiplierCPU(map, normalSpeed, localArea, Dt, c11/c22); 
	localTension = map.WENO5RK3Extend(localTension, 100);
	residual = map.WENO5RK3Extend(residual, 100);

	% calculate tension and pressure
	normalSpeed = NormalSpeedBend - map.MeanCurvature .* localTension;
	%s1 = - map.surfaceIntegral(map.MeanCurvature .* NormalSpeedBend) ...
	s1 =  map.surfaceIntegral(map.MeanCurvature .* normalSpeed ) ...
		 - (InitialArea - CurrentArea) / Dt;
	s2 = map.surfaceIntegral(normalSpeed) + (InitialVolume - CurrentVolume) / Dt;
	
	TP = [c11,c12;c21,c22]\[s1;s2];
	Tension = TP(1);
	Pressure = TP(2);

	% now calculate normal Speed that preserves total area and volume
	normalSpeed = Tension * map.MeanCurvature - normalSpeed + Pressure;

	% time step c field with first order Euler scheme
	TotalC = map.surfaceIntegral(localArea);
	%DiffC = 100 * (TotalC - InitialArea) / InitialArea;
	DiffC = 100 * (TotalC - CurrentArea) / InitialArea;
	residual = residual - map.MeanCurvature .* (map.MeanCurvature .* Tension + Pressure);
	%localArea = localArea .* ( 1.0 + Dt * map.MeanCurvature .* normalSpeed);
	localArea = timeStepLocalArea(map, localArea, localTension, residual, Dt, 0.5);
	localArea = map.WENO5RK3Extend(localArea, 100);

	% time step level set function
	% addtion of localTnesion causes total area and volume to change
	normalSpeedSmoothed = smoothGMRES(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);
	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance

	% extend local Area field away from the new surface position
	localArea = map.WENO5RK3Extend(localArea, 100);

	% %5d is important to avoid mess in ssh matlab sesssion 
	fprintf('iter: %5d, ene: %4.5f, ar: %4.5f, vol: %4.5f, rd: %4.5f, c: %4.5f\n', ...
			i, c11, DiffArea, DiffVolume, ReduceVolume, DiffC)
	%fprintf('iter: %5d, ene: %4.5f, DiffC: %4.5f\n', i, c11, DiffC)

	if mod(i,1)==pi 
		timeStr = [sprintf('%04d: %0.5e, %0.5f', i,time,c11)];

		clf(FIG)

		%map.plotSurface(0,1,'Green','black');textZ = gather(map.GD3.zmin);
		%map.plotField(0,normalSpeedSmoothed,0.5)
		map.plotField(0,localArea,0.5)
		%map.plotField(0,localTension,0.5)
		ax = gca;
		ax.Visible = 'off';
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
	end


end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt + Alpha * Dt * map.GD3.LBiLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

	[normalSpeedSmoothed,~,~,~,~] = gmres(Op, S, [], 1e-12, 300, @mfun);
	%normalSpeedSmoothed = gmres(Op, S, [], 1e-6, 300, @mfun);
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

% calculate local Lagrange multiplier
function [localTension,residual] = localLagrangeMultiplier(map, normalSpeed, localArea, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Lxx + map.GD3.Lyy + map.GD3.Lzz - ...
			(   map.GD3.SparseDiag(map.Nx .^2) * map.GD3.Lxx + ...
				map.GD3.SparseDiag(map.Ny .^2) * map.GD3.Lyy + ...
				map.GD3.SparseDiag(map.Nz .^2) * map.GD3.Lzz + ...
				map.GD3.SparseDiag(map.Nx .* map.Ny) * map.GD3.Lxy * 2 + ...
				map.GD3.SparseDiag(map.Ny .* map.Nz) * map.GD3.Lyz * 2 + ...
				map.GD3.SparseDiag(map.Nz .* map.Nx) * map.GD3.Lzx * 2 ...
			) ...
		- map.GD3.SparseDiag(map.MeanCurvature.^2) ;
	% right hand side
	rhs = map.MeanCurvature .* normalSpeed + ((localArea - 1) ./ localArea) / Dt;
	S = reshape(rhs, [map.GD3.NumElt, 1]);

	% it seems that restart about {10,...,15} gives best speed
	[localTension,~,~,~,~] = gmres(Op, S, 11, 1e-6, 300, @mfun);
	residual = Op * localTension - S;

	%localTension = gmres(Op, S, 10, 1e-13, 300, @mfun);
	localTension = reshape(localTension, map.GD3.Size);
	residual = reshape(residual, map.GD3.Size);

	% preconditioner, Alpha = c11/c22 as the mean squared MeanCurvature
	function y = mfun(S)
		fftS = fftn(reshape(S,map.GD3.Size));
		fftS = - fftS ./ (Alpha + map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2);
		y = real(ifftn(fftS));
		y = reshape(y, [map.GD3.NumElt, 1]);
	end

end

% this version uses ilu decomposion to precondition gmeres
% while on smaller grid [64,64,32], the speed is comparable to the GPU version
% on largerer grid [128,128,64] the speed is 30 times slower
% GPU code scales much much better than cpu code
function localTension = localLagrangeMultiplierCPU(map, normalSpeed, localArea, Dt, Alpha)

	tic
	% operator to be solved
	Op = map.GD3.Lxx + map.GD3.Lyy + map.GD3.Lzz - ...
			(   map.GD3.SparseDiag(map.Nx .^2) * map.GD3.Lxx + ...
				map.GD3.SparseDiag(map.Ny .^2) * map.GD3.Lyy + ...
				map.GD3.SparseDiag(map.Nz .^2) * map.GD3.Lzz + ...
				map.GD3.SparseDiag(map.Nx .* map.Ny) * map.GD3.Lxy * 2 + ...
				map.GD3.SparseDiag(map.Ny .* map.Nz) * map.GD3.Lyz * 2 + ...
				map.GD3.SparseDiag(map.Nz .* map.Nx) * map.GD3.Lzx * 2 ...
			) ...
		- map.GD3.SparseDiag(map.MeanCurvature.^2) ;
	% right hand side
	rhs = map.MeanCurvature .* normalSpeed + ((localArea - 1) ./ localArea) / Dt;
	S = reshape(rhs, [map.GD3.NumElt, 1]);

	OpCPU = gather(Op);
	SCPU = gather(S);

	[L,U] = ilu(OpCPU, struct('type','nofill','milu','row'));
	localTensionCPU = gmres(OpCPU, SCPU, 10, 1e-6, 300, L, U);

	localTension = reshape(gpuArray(localTensionCPU), map.GD3.Size);

	toc

end

% time step c field
function NewC = timeStepLocalArea(map, localArea, localTension, residual, Dt, Alpha)

	Advection = zeros(map.GD3.Size, 'gpuArray');
	[vx,vy,vz] = map.GD3.Gradient(localTension);
	Advection = feval(map.advection_step,Advection,vx,vy,vz,localArea,...
			map.GD3.mrows,map.GD3.ncols,map.GD3.lshts,...
			map.GD3.Dx,map.GD3.Dy,map.GD3.Dz);

	rhs = 1.0 - Dt*(Alpha * map.GD3.Laplacian(localArea) + Advection + residual.* localArea);

	% now use fft to solve for newC
	fftS = fftn(rhs);
	fftS = fftS ./ (1.0 + Alpha * Dt * ...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2) );
	NewC = real(ifftn(fftS));


end












