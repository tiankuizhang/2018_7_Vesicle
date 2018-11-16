% test order of convergence of the new surface redistance scheme

for ii = [64,96,128,196]
	test_surface_redistance(ii);
end

function test_surface_redistance(N)
	r0 = 0.6;
	xmin = -1.0; xmax = 1.0;
	
	xv = linspace(xmin,xmax,N);
	yv = xv;
	zv = xv;
	
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	grid = SD.GD3(x,y,z);
	[azimuth,elevation,r] = cart2sph(x,y,z);
	
	F= sqrt(x.^2 + y.^2 + z.^2) - r0;
	map = SD.SDF3(grid,x,y,z,F);
	map.setDistance
	map.setCalculusToolBox4;
	
	map.A = r0 .* sin(elevation);
	
	G_exact = r0 .* elevation;
	%G_nu = map.WENO5RK3ClosetPointSurfaceRedistance(map.A, 300, 20);
	%G_nu = map.WENORK3ClosetPointSurfaceRedistance(map.A, 300, 20);
	G_nu = map.ENORK2ClosetPointSurfaceRedistance(map.A, 300, 20);
	map.A = G_nu;
	map.AsetCalculusToolBox4;
	
	GeodesicCurvature_exact = tan(elevation) ./ r0;
	mask = abs(GeodesicCurvature_exact) > 4.0;
	
	GeodesicCurvature_exact(mask) = 0.0;
	map.GeodesicCurvature(mask) = 0.0;
	diff = GeodesicCurvature_exact - map.GeodesicCurvature;
	
	diff = map.WENO5RK3Extend(G_exact-G_nu,100);
	
	E = sqrt(map.surfaceIntegral(diff.^2));
	fprintf('N: %5d, E: %5.3e \n',N, E)
end

