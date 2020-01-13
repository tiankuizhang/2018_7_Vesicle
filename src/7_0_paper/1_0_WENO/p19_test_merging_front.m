% test the merging sphere example from reviewer2


%	Nx = 16; Ny = 16; Nz = 24;
%	Nx = 32; Ny = 32; Nz = 48;
	Nx = 64; Ny = 64; Nz = 96;
%	Nx = 128; Ny = 128; Nz = 192;

%	Nx = 16; Ny = 16; Nz = 32;
%	Nx = 32; Ny = 32; Nz = 64;
%	Nx = 64; Ny = 64; Nz = 128;
%	Nx = 128; Ny = 128; Nz = 256;

	r = .6;
	xmin = -1.0; xmax = 1.0;
	xv = linspace(xmin, xmax, Nx);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	%yv = xv; zv = xv;

	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	grid = SD.GD3(x,y,z);

	z0 = 0.5 * r;
	F1 = sqrt( x.^2 + y.^2 + (z-z0).^2 ) - r; 
	F2 = sqrt( x.^2 + y.^2 + (z+z0).^2 ) - r; 
	F_exact = min(F1,F2);
	%F_exact = F1;

	f1 = x.^2 + y.^2 + (z-z0).^2 - r^2;
	f2 = x.^2 + y.^2 + (z+z0).^2 - r^2;

	%F = f1; F(maskF2) = f2(maskF2);
	F = min(f1,f2);

	map = SD.SDF3(grid,x,y,z,F);

%	figure
%	subplot(1,2,1)
%	map.plot()

	map.setDistance
	map.F = map.WENO5RK3Reinitialization(map.F,1000);
%	map.F = map.WENORK3Reinitialization(map.F,1000);
%	map.F = map.ENORK2Reinitialization(map.F,1000);
%	subplot(1,2,2)
%	map.plot()

	mask = abs(map.F) < 1.5*map.GD3.Ds;
%	map.setCalculusToolBox4;
	map.setCalculusToolBox;

	% error in level set function
	DF = map.F - F_exact;
	DF = map.WENORK3Extend(DF,100);
	EF1 = map.surfaceIntegral(abs(DF));
	EF2 = sqrt(map.surfaceIntegral(DF.^2));
	EFM = max(abs(DF(mask)));

%	map.plotField(0,DF,0)

	fprintf('%03d:    F1: %5.3e,    F2: %5.3e,    FM: %5.3e , z0:%.2f \n',Nx,EF1,EF2,EFM,z0)






























