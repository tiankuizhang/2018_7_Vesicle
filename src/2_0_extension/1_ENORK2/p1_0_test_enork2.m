% test ENO-RK2 reinitialization scheme

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
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-0.6;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

obj = map;

Fgpu = obj.F;

C = obj.GD3.Z;

xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
		xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

% calculate the extend velocity upwindly
fx = zeros(obj.GD3.Size, 'gpuArray');
fy = zeros(obj.GD3.Size, 'gpuArray');
fz = zeros(obj.GD3.Size, 'gpuArray');

[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
		fx, fy, fz, Fgpu, xpr, xpl, ypf, ypb, zpu, zpd, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

fgradient = sqrt(fx.^2+fy.^2+fz.^2);
flat = fgradient < 1e-6;
fgradient(flat) = 1e-6;

nx = fx ./ fgradient;
ny = fy ./ fgradient;
nz = fz ./ fgradient;

vx = sign(Fgpu) .* nx;
vy = sign(Fgpu) .* ny;
vz = sign(Fgpu) .* nz;

% boundary nodes
Boundary = ( obj.F .* obj.F(obj.GD3.oxo) < 0 ) | ...
		   ( obj.F .* obj.F(obj.GD3.oXo) < 0 ) | ...
		   ( obj.F .* obj.F(obj.GD3.yoo) < 0 ) | ...
		   ( obj.F .* obj.F(obj.GD3.Yoo) < 0 ) | ...
		   ( obj.F .* obj.F(obj.GD3.ooz) < 0 ) | ...
		   ( obj.F .* obj.F(obj.GD3.ooZ) < 0 ) ; 

% closest points on boundary from boundary nodes 
x_shift = obj.GD3.X(Boundary) - nx(Boundary) .* obj.F(Boundary);
y_shift = obj.GD3.Y(Boundary) - ny(Boundary) .* obj.F(Boundary);
z_shift = obj.GD3.Z(Boundary) - nz(Boundary) .* obj.F(Boundary);

C(Boundary) = interp3(obj.GD3.X, obj.GD3.Y, obj.GD3.Z, C, ...
		x_shift, y_shift, z_shift, 'linear');


step = zeros(obj.GD3.Size, 'gpuArray');
deltat = 0.3 * min(obj.GD3.Dx, min(obj.GD3.Dy,obj.GD3.Dz));

toc;

tic;
for i=1:30
	step = feval(obj.ENORK2_extend_step, step, C, vx, vy, vz, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	step(Boundary) = 0;
	Ctmp = C - deltat * step;
	step = feval(obj.ENORK2_extend_step, step, Ctmp, vx, vy, vz, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	step(Boundary) = 0;
	C = (C + Ctmp - deltat * step) / 2;
end


toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
