function test_error(N)
	% create a 3D grid
	xv = linspace(-1,1,N);
	yv = xv;
	zv = xv;
	
	[x,y,z] = meshgrid(xv,yv,zv);
	
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	
	grid = SD.GD3(x,y,z);
	
	% create a SDF3 instance
	Radius = 0.6;
	fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-Radius;
	
	F = fun(x, y, z);
	
	map = SD.SDF3(grid, x, y, z, F);
	
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% test convergence and accuracy of the smeared Dirac_Delta function 
	% and Heaviside function
	% following methods from towers
	
	tic
	Primal = max(map.F,0);
	
	Px = map.GD3.Fx(Primal);
	Py = map.GD3.Fy(Primal);
	Pz = map.GD3.Fz(Primal);
	
	P_lap = map.GD3.Laplacian(Primal);
	
	Fx = map.GD3.Fx(map.F);
	Fy = map.GD3.Fy(map.F);
	Fz = map.GD3.Fz(map.F);
	F_lap = map.GD3.Laplacian(map.F);
	
	mag_grad = sqrt(Fx.^2+Fy.^2+Fz.^2);
	
	dot_Pirmal_F = Px.*Fx + Py.*Fy + Pz.*Fz;
	
	Heaviside = dot_Pirmal_F ./ mag_grad.^2;
	
	Dirac_Delta = P_lap ./ mag_grad.^2 - dot_Pirmal_F .* F_lap ./ mag_grad.^4;
	toc
	
	mask1 = abs(map.F) < 2*map.GD3.Dx;
	mask2 = map.F < 2*map.GD3.Dx;
	
	NuArea = sum(Dirac_Delta(mask1)) * map.GD3.Dx.^3;
	ThArea = 4*pi*Radius.^2;
	AreaRelativeError = abs((ThArea-NuArea)/ThArea)
	
	NuVol = sum(1-Heaviside(mask2)) * map.GD3.Dx.^3;
	ThVol = 4*pi*Radius.^3/3;
	VolRelativeError = abs((ThVol-NuVol)/ThVol)

end	



 












