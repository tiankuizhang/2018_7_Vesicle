% solve surface diffusion problem with an semi-implicit scheme
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BiLaplacian = map.GD3.LBiLaplacian;
%Opgpu = 0.001 * BiLaplacian + map.GD3.Idt;
%Opcpu = gather(Opgpu);
%Scpu = rand(128^3, 1);
%%Scpu = rand(64^3, 1);
%Sgpu = gpuArray(Scpu);
%
%[L, U] = ilu(Opcpu, struct('type','nofill','milu','row'));
%LUgpu = gpuArray(L*U);
%
%
%[L, U] = ilu(B, struct('type','nofill','milu','row'))
% F_new = gmres(B, S,[], 1e-12, 300, L,U)
% F_new = bicgstab(B, S, 1e-12, 300, L,U)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% name and size of figure
FIG = figure('Name','Surface Diffusion of perturbed cylinder','Position',[10 10 1600 800])

time = 0;

MaxResolvedCurvature = 2.0 / map.GD3.Ds;

textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

% CFL number should be small enough to avoid oscillation
% 0.1 is a safe choice for forward Euler time step
CFLNumber = 1;

for i=1:20000
	
	map.GPUsetCalculusToolBox
	mask = abs(map.F)<2*map.GD3.Ds;
	
	map.MeanCurvature = sign(map.MeanCurvature) .* ...
		min(MaxResolvedCurvature, abs(map.MeanCurvature));
	MaxCurvatureBeforeExtend = max(abs(map.MeanCurvature(mask)));

	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	%MeanCurvature = map.ENORK2Extend(map.MeanCurvature,100);
	MaxCurvature = max(abs(MeanCurvature(mask)));

	% surface Laplacian of mean curvature
	
	MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian(MeanCurvature); 
	NormalSpeed = - map.WENORK3Extend(MeanCurvatureSurfaceLaplacian .* map.FGradMag,100);

	MaxSpeed = max(abs(NormalSpeed(mask)));
	
	Dt = CFLNumber * map.GD3.Ds / MaxSpeed;
	time = time + Dt;
	
	fprintf('iter: %d, time: %4.10f, MaxCrv: %4.5f, MaxCrv: %4.5f Dt: %4.5f\n', i,time, ...
			MaxCurvatureBeforeExtend/MaxResolvedCurvature, ...
			MaxCurvature/MaxResolvedCurvature,Dt/map.GD3.Ds^4)

	% now use FFT to smooth NormalSpeed
	normalSpeedFFT = fftn(NormalSpeed);
	normalSpeedFFTsmoothed = normalSpeedFFT ./ ( 1 + 0.5 * Dt *	...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
	normalSpeedSmoothed = real(ifftn(normalSpeedFFTsmoothed));

	normalSpeedSmoothed = map.WENORK3Extend(normalSpeedSmoothed, 100);


	%map.F = map.F - Dt * NormalSpeed;
	map.F = map.F - Dt * normalSpeedSmoothed;

	if mod(i,1)==0
		timeStr = [sprintf('%04d: %0.5e', i,time)];

		clf(FIG)

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

	end

	% the system should have enough time for scission between two reinitialization step
	% otherwise, scission will never happen
	% it is best to avoid reinitialization all together
	if mod(i,50)==0
		%map.F = map.WENORK3Reinitialization(map.F, 100);
	end

end

