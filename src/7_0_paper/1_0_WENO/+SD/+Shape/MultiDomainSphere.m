% initialize a sphere with multiple domains
% domain = [az0,ele0,alpha; ...; ...]

function [x,y,z,F,A,volume] = MultiDomainSphere(boundingBox, Size, radius, ...
		reducedVolume, domain)

	% create a meshgrid
	Nx = Size(1); Ny = Size(2); Nz = Size(3);
	xmin = boundingBox(1); xmax = boundingBox(2);
	xv = linspace(xmin,xmax,Nx); dx = xv(2) - xv(1);
	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x); y = gpuArray(y); z = gpuArray(z);

	r = sqrt(x.^2+y.^2+z.^2);
	F = r - radius;
	volume = reducedVolume * (4*pi/3) * radius^3;

	% construct A function
	rx = x./r; ry = y./r; rz = z./r; % unit radial vector
	[n,~] = size(domain);
	A = Inf * ones(Size);
	for i=1:n
		az0 = domain(i,1); ele0 = domain(i,2); alpha = domain(i,3);
		A = min(A,SphericalCap(rx,ry,rz,az0,ele0,alpha));
	end

end


function A = SphericalCap(rx,ry,rz,az0,ele0,alpha)

	% unit vector in (az0,ele0) direction
	rx0 = cos(ele0) * cos(az0);
	ry0 = cos(ele0) * sin(az0);
	rz0 = sin(ele0);

	A = cos(alpha) - (rx.*rx0 + ry.*ry0 + rz.*rz0) ;


end
