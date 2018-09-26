% return an Ellipsoid with the specified reduced volume
% the default half axis are a=1,b=1,c from interpolant data

function [x,y,z,f] = Ellipsoid(N, reducedVolume)

	% create a meshgrid
	xmin = -1.0;
	xmax = 1.0;
	xcenter = (xmin + xmax) / 2;
	xv = linspace(xmin, xmax, N);
	yv = xv;
	zv = xv;

	[x,y,z] = meshgrid(xv,yv,zv);

	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);

	% load reducedvolume vs c data if exists, otherwise calculate it
	FILE = 'EllipsoidReducedVolume.mat';
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

	% create a level set function whose zero contour is the desired shape
	a = 0.35 * (xmax - xmin);
	b = a;
	c = a * cq;

	f = sqrt((x-xcenter).^2 ./ a^2 + (y-xcenter).^2 ./ b^2 + (z-xcenter).^2 ./ c^2) - 1;

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
