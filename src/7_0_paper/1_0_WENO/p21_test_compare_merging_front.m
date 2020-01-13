
r = 0.6;
N = 19;
NxList = [16,32,64,128];
%NxList = [128];
raList = linspace(0.5,1.5,N);

EM = [];
Dist = 2*r*(raList - 1);

for j = 1:4
	tmp = [];
	for i = 1:N
		[EF1,EF2,EFM] = CalculateError(NxList(j),r,raList(i));
		tmp = [tmp EFM];
	end
	EM = [EM; tmp];
end

EM = gather(EM);
save('mergingSphere','EM','Dist')

function [EF1,EF2,EFM] = CalculateError(N,r,ra)
	Nx = N; Ny = N; Nz = 2*N;

	r = .6;
	xmin = -1.0; xmax = 1.0;
	xv = linspace(xmin,xmax,N);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	grid = SD.GD3(x,y,z);
	
	z0 = ra * r;
	F1 = sqrt( x.^2 + y.^2 + (z-z0).^2 ) - r; 
	F2 = sqrt( x.^2 + y.^2 + (z+z0).^2 ) - r; 
	F_exact = min(F1,F2);


	f1 = x.^2 + y.^2 + (z-z0).^2 - r^2;
	f2 = x.^2 + y.^2 + (z+z0).^2 - r^2;
	F = min(f1,f2);
	
	map = SD.SDF3(grid,x,y,z,F);
	map.setDistance
	map.F = map.WENO5RK3Reinitialization(map.F,1000);
	%map.F = map.WENORK3Reinitialization(map.F,1000);
	map.setCalculusToolBox4;

	mask = abs(map.F) < 1.5*map.GD3.Ds;
	
	% error in level set function
	DF = map.F - F_exact;
	DF = map.WENORK3Extend(DF,100);
	EF1 = map.surfaceIntegral(abs(DF));
	EF2 = sqrt(map.surfaceIntegral(DF.^2));
	EFM = max(abs(DF(mask)));
	
	fprintf('%03d: F1: %5.3e, F2: %5.3e, FM: %5.3e, Di: %.3f \n',N,EF1,EF2,EFM,2*z0-2*r)
end
