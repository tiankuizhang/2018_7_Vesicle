% inputs:	boundingBox -- [xmin, xmax]
%			size -- [Nx, Ny, Nz]
%			radius, reduced volume, Lipid disordered phase area ratio
% outputs:	 meshgrid --- x,y,z
%			level set function and auxilary level set function -- F,A
%			expected volume -- v

function [x,y,z,F,A,volume] = BiphaseSphere(boundingBox, Size, radius, reducedVolume, raLd)

	% create a meshgrid
	Nx = Size(1); Ny = Size(2); Nz = Size(3);
	xmin = boundingBox(1); xmax = boundingBox(2);
	xv = linspace(xmin,xmax,Nx); dx = xv(2) - xv(1);
	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x); y = gpuArray(y); z = gpuArray(z);

	F = sqrt(x.^2+y.^2+z.^2) - radius;
	
	z0 = radius * (1 - 2*raLd); % height of phase boundary, z>z0 represents Ld phase
	if z0 ~= 0 
		A = z0 - sign(z0) * z;
	else
		A = - z;
	end

	volume = reducedVolume * (4*pi/3) * radius^3;

end
