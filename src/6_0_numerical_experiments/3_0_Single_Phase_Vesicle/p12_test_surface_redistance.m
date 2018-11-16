% test order of convergence of the new surface redistance scheme

%for ii = [64,96,128]
for ii = [64,96,128]
	test_surface_redistance(ii);
end

%N = 64;
%N = 128;
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
	%G_nu = map.WENO5RK3ClosetPointSurfaceRedistance(map.A, 100, 100);
	%G_nu = map.WENORK3ClosetPointSurfaceRedistance(map.A, 100, 100);
	G_nu = map.ENORK2ClosetPointSurfaceRedistance(map.A, 100, 100);
	map.A = G_nu;
	map.AsetCalculusToolBox4;
	
	GeodesicCurvature_exact = tan(elevation) ./ r0;

	mask = abs(GeodesicCurvature_exact) > 2;
	%GeodesicCurvature_exact(mask) = 0.0;
	%map.GeodesicCurvature(mask) = 0.0;

	diff_G = G_nu - G_exact;
	diff_G = map.WENO5RK3Extend(diff_G,100);
	E_G = sqrt(map.LineIntegral(diff_G.^2));
	diff_G(mask) = 0.0;
	E_G_s = sqrt(map.surfaceIntegral(diff_G.^2));

	diff_c = GeodesicCurvature_exact - map.GeodesicCurvature;
	diff_c = map.WENO5RK3Extend(diff_c,100);
	E_c = sqrt(map.LineIntegral(diff_c.^2));
	diff_c(mask) = 0.0;
	E_c_s = sqrt(map.surfaceIntegral(diff_c.^2));
	
	fprintf('N: %5d, E_G: %5.3e, E_c: %5.3e, E_G_s: %5.3e, E_c_s: %5.3e \n',N, E_G, E_c,E_G_s,E_c_s)
end

