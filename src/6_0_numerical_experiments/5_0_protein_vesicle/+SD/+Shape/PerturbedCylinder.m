% return a perturbed cylinder with a specified beta value between 0.6 and 0.9

function [x,y,z,f] = PerturbedCylinder(Size, Beta)

	% create a meshgrid
	Nx = Size(1);
	Ny = Size(2);
	Nz = Size(3);

	D = 0.5;
	Alpha = 2;

	xmin = -2*pi - D;
	xmax =  2*pi + D;
	xv = linspace(xmin, xmax, Nx) / Alpha;
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	zmin = zv(1);
	zmax = zv(end);

	[x,y,z] = meshgrid(xv,yv,zv);

	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);

	% create a perturbed cylinder
	C1 = 0.5 * pi / Alpha; % maximum radius
	C2 = Beta * C1 /2; % C1 - 2*C2 will be the minmum radius

	F1 = sqrt(x.^2+y.^2) - (C1-C2*(cos(Alpha * z) + 1));
	F2 = max(z-pi/Alpha,-z-pi/Alpha);

	f = max(F1,F2);

	F3 = sqrt(x.^2+y.^2+(z-pi/Alpha).^2) - C1;
	F4 = sqrt(x.^2+y.^2+(z+pi/Alpha).^2) - C1;

	F5 = min(F3,F4);

	f = min(f,F5);

end
