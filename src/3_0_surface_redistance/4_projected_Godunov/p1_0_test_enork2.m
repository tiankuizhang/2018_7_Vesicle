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

%AF = obj.ENORK2Extend(obj.GD3.Z,100);
AF = obj.GD3.Z-0.3;
OldAF = AF;

xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

% AF instead of Fgpu should be used!!!!
[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
		xpr, xpl, ypf, ypb, zpu, zpd, AF, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

% calculate the normal to the surface upwindly
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

% construct a local coordinate at each point
ax = zeros(obj.GD3.Size, 'gpuArray');
ay = zeros(obj.GD3.Size, 'gpuArray');
az = zeros(obj.GD3.Size, 'gpuArray');

bx = zeros(obj.GD3.Size, 'gpuArray');
by = zeros(obj.GD3.Size, 'gpuArray');
bz = zeros(obj.GD3.Size, 'gpuArray');

[ax,ay,az,bx,by,bz] = feval(obj.surface_coordinate,ax,ay,az,bx,by,bz,nx,ny,nz,...
		obj.GD3.mrows,obj.GD3.ncols,obj.GD3.lshts);

% calculate the sign of the original auxiliary level set function
Sign = zeros(obj.GD3.Size, 'gpuArray');
Sign(AF>0) = 1.;
Sign(AF<0) = -1.;

%
minx = min(xpr,xpl);
miny = min(ypf,ypb);
minz = min(zpu,zpd);
deltat = 0.3 * min(minx, min(miny,minz));

step = zeros(obj.GD3.Size, 'gpuArray');

tic
for i=1:50
	step = feval(obj.ENORK2_surface_redistance_step,step,AF,Sign,deltat,...
		   		nx,ny,nz,ax,ay,az,bx,by,bz,xpr,xpl,ypf,ypb,zpu,zpd,...	
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	AFtmp = AF - step;
	step = feval(obj.ENORK2_surface_redistance_step,step,AFtmp,Sign,deltat,...
		   		nx,ny,nz,ax,ay,az,bx,by,bz,xpr,xpl,ypf,ypb,zpu,zpd,...	
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	AF = (AF + AFtmp - step) / 2;

end
toc

AF = obj.ENORK2Extend(AF,100);

%	step = feval(obj.ENORK2_surface_redistance_step,step,AF,Sign,deltat,nx,ny,nz,...
%		   		xpr,xpl,ypf,ypb,zpu,zpd,...	
%		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
%				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure

subplot(1,2,1)
obj.plotIsoField([0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7],OldAF)

subplot(1,2,2)
obj.plotIsoField([0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7],AF)




























