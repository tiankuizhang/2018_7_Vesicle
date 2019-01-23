% return an Ellipsoid with the specified reduced volume
% the default half axis are a=1,b=1,c for oblate
% and a=1, b=c for prolate 

function [x,y,z,f,a,b,c] = Ellipsoid(Size, reducedVolume, TYPE)

	ratio = 0.4;
	% create a meshgrid
	Nx = Size(1);
	Ny = Size(2);
	Nz = Size(3);

	xmin = -1.0;
	xmax = 1.0;
	xv = linspace(xmin, xmax, Nx);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	zmin = zv(1);
	zmax = zv(end);

	[x,y,z] = meshgrid(xv,yv,zv);

	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);

	% load reducedvolume vs c data if exists, otherwise calculate it
	if TYPE=='Oblate' || TYPE=='O' || TYPE=='o'
		FILE = fullfile('+SD','+Shape','EllipsoidOblateReducedVolume.mat');
		if exist(FILE)
			load(FILE)
		else
			c = 0.01:0.01:0.99;
			rv = zeros(size(c));
			for ii = 1:length(c);
				rv(ii) = ReducedVolume(1,1,c(ii));
			end
			save(FILE,'rv','c')
		end
		cq = interp1(rv, c, reducedVolume); 
		a = ratio * (xmax - xmin);
		b = a;
		c = a * cq;
		f = sqrt(x.^2 ./ a^2 + y.^2 ./ b^2 + z.^2 ./ c^2) - 1;
	end

	if TYPE=='Prolate' || TYPE=='P' || TYPE=='p'
		FILE = fullfile('+SD','+Shape','EllipsoidProlateReducedVolume.mat');
		if exist(FILE)
			load(FILE)
		else
			c = 0.01:0.01:0.99;
			rv = zeros(size(c));
			for ii = 1:length(c);
				rv(ii) = ReducedVolume(1,c(ii),c(ii));
			end
			save(FILE,'rv','c')
		end
		cq = interp1(rv, c, reducedVolume); 
		a = ratio * (zmax - zmin);
		b = a * cq;
		c = b;
		f = sqrt(x.^2 ./ b^2 + y.^2 ./ c^2 + z.^2 ./ a^2) - 1;
	end



end

% calculate reduced volume for an ellipsoid with half axis a>=b>c 
function rv = ReducedVolume(a,b,c)
	Volume = 4 * pi * a * b * c / 3;
	
	phi = atan2(sqrt(a^2-c^2),c);
	cos_phi = c/a;
	sin_phi = sqrt(a^2-c^2)/a;
	k = a * sqrt(b^2-c^2) / b / sqrt(a^2-c^2);
	
	Area = 2*pi*c^2 + 2*pi*a*b*(ellipticE(phi,k)*sin_phi^2 + ellipticF(phi,k)*cos_phi^2) / sin_phi;
	
	rv = (3*Volume/4/pi) * (4*pi/Area)^(3/2);
end
