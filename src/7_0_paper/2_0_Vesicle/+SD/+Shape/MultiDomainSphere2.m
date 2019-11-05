% initialize a sphere with multiple domains
% domain = [az0,ele0,alpha,beta; ...; ...]

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
		az0 = domain(i,1); ele0 = domain(i,2); 
		alpha = domain(i,3); beta = domain(i,4);
		[A0,F0] = SphericalCap(F,radius,x,y,z,rx,ry,rz,az0,ele0,alpha,beta);
		A = min(A,A0);
		F = max(F,F0);
	end

end


function [A0,F0] = SphericalCap(F,radius,x,y,z,rx,ry,rz,az0,ele0,alpha,beta)

	% unit vector in (az0,ele0) direction
	rx0 = cos(ele0) * cos(az0);
	ry0 = cos(ele0) * sin(az0);
	rz0 = sin(ele0);

	A0 = cos(alpha) - (rx.*rx0 + ry.*ry0 + rz.*rz0) ;

	if beta > 0
		r0 = radius * sin(alpha) / sin(beta); % radius of F0
		d = radius * cos(alpha) + r0 * cos(beta); % distance between two sphere
		cx0 = rx0 * d; % center for F0
		cy0 = ry0 * d;
		cz0 = rz0 * d;
	
		F0 = r0 - sqrt( (x-cx0).^2 + (y-cy0).^2 + (z-cz0).^2 ); 
	else
		F0 = F;
	end

end
