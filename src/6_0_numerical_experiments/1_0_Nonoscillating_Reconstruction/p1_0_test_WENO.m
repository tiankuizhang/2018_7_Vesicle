% create a 3D grid
xv = linspace(-1,1,64);
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
map.A = z;

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct WENO derivatives of Field 

obj = map;
Field = x.^2 + y.^3 + z.^4 + x.*y.*z;
% Field = z;

% WENO derivatives in x direction 
DFx = (Field - Field(obj.GD3.oxo)) / obj.GD3.Dx;

v1 = circshift(DFx, [ 0, 2, 0]);
v2 = circshift(DFx, [ 0, 1, 0]);
v3 = DFx;
v4 = circshift(DFx, [ 0,-1, 0]);
v5 = circshift(DFx, [ 0,-2, 0]);
v6 = circshift(DFx, [ 0,-3, 0]);

WENO_back_x = WENO_Derivative(v1, v2, v3, v4, v5);
WENO_fore_x = WENO_Derivative(v6, v5, v4, v3, v2);

% WENO derivatives in y direction 
DFy = (Field - Field(obj.GD3.yoo)) / obj.GD3.Dy;

v1 = circshift(DFy, [ 2, 0, 0]);
v2 = circshift(DFy, [ 1, 0, 0]);
v3 = DFy;
v4 = circshift(DFy, [-1, 0, 0]);
v5 = circshift(DFy, [-2, 0, 0]);
v6 = circshift(DFy, [-3, 0, 0]);

WENO_back_y = WENO_Derivative(v1, v2, v3, v4, v5);
WENO_fore_y = WENO_Derivative(v6, v5, v4, v3, v2);

% WENO derivatives in x direction 
DFz = (Field - Field(obj.GD3.ooz)) / obj.GD3.Dz;

v1 = circshift(DFz, [ 0, 0, 2]);
v2 = circshift(DFz, [ 0, 0, 1]);
v3 = DFz;
v4 = circshift(DFz, [ 0, 0,-1]);
v5 = circshift(DFz, [ 0, 0,-2]);
v6 = circshift(DFz, [ 0, 0,-3]);

WENO_back_z = WENO_Derivative(v1, v2, v3, v4, v5);
WENO_fore_z = WENO_Derivative(v6, v5, v4, v3, v2);

N = 100;
tic
for i=1:N
	[xl,xr,yb,yf,zd,zu] = map.GD3.GPUWENODerivative(Field);
end
toc

tic
for i=1:N
	[xl,xr,yb,yf,zd,zu] = map.GD3.WENODerivative(Field);
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WENO_D = WENO_Derivative(v1, v2, v3, v4, v5)

	c11 = 1./3.;	c12 = -7./6.;	c13 = 11./6.;
	c21 = -1./6.;	c22 = 5./6.;	c23 = 1./3.;
	c31 = 1./3;		c32 = 5./6.;	c33 = -1./6.;

	% different choices of ENO derivatives
	phi1 = c11 * v1 + c12 * v2 + c13 * v3;
	phi2 = c21 * v2 + c22 * v3 + c23 * v4;
	phi3 = c31 * v3 + c32 * v4 + c33 * v5;

	% smoothness parameter
	S1 = 13./12.*(v1 - 2*v2 + v3).^2 + 1./4.*(v1 - 4*v2 + 3*v3).^2;
	S2 = 13./12.*(v2 - 2*v3 + v4).^2 + 1./4.*(v2 - v4).^2;
	S3 = 13./12.*(v3 - 2*v4 + v5).^2 + 1./4.*(3*v3 - 4*v4 + v5).^2;
	
	epsilon = 1e-6;
	alpha1 = 0.1 ./ (S1 + epsilon).^2;
	alpha2 = 0.6 ./ (S2 + epsilon).^2;
	alpha3 = 0.3 ./ (S3 + epsilon).^2;
	
	% weights for each choice
	omega1 = alpha1 ./ (alpha1 + alpha2 + alpha3);
	omega2 = alpha2 ./ (alpha1 + alpha2 + alpha3);
	omega3 = alpha3 ./ (alpha1 + alpha2 + alpha3);
	
	WENO_D = omega1 .* phi1 + omega2 .* phi2 + omega3 .* phi3;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






























