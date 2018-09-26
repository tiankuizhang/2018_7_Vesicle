% find parameters for ellipsoid of different reduced volume

a = 1;
b = 1;
c = 0.01:0.01:0.99;

rv = zeros(size(c));
tic
for ii = 1:length(c)
	rv(ii) = ReducedVolume(a,b,c(ii));
end
toc

save('ReducedVolume.mat','rv','c')

tic
load('ReducedVolume.mat')
toc

cq = interp1(rv, c, 0.5)

ReducedVolume(a,b,cq)


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
