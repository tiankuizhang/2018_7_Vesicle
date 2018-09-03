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
% solve the surface redistancing equation with the closet point method

obj = map;
original = map.GD3.Z + 0.3;
oldC = obj.ENORK2Extend(original, 100);
C = oldC;

% calculate distance to the neighborhood node without crossing the interface
xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
		xpr, xpl, ypf, ypb, zpu, zpd, oldC, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

minx = min(xpr,xpl);
miny = min(ypf,ypb);
minz = min(zpu,zpd);
deltat = 0.3 * min(minx, min(miny,minz));

% calculate the extend velocity upwindly
fx = zeros(obj.GD3.Size, 'gpuArray');
fy = zeros(obj.GD3.Size, 'gpuArray');
fz = zeros(obj.GD3.Size, 'gpuArray');

[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
		fx, fy, fz, obj.F, xpr, xpl, ypf, ypb, zpu, zpd, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

fgradient = sqrt(fx.^2+fy.^2+fz.^2);
flat = fgradient < 1e-6;
fgradient(flat) = 1e-6;

nx = fx ./ fgradient;
ny = fy ./ fgradient;
nz = fz ./ fgradient;

band = abs(map.F) < 6 * map.GD3.Dx;

% closest points on boundary from boundary nodes 
x_shift = obj.GD3.X(band) - nx(band) .* obj.F(band);
y_shift = obj.GD3.Y(band) - ny(band) .* obj.F(band);
z_shift = obj.GD3.Z(band) - nz(band) .* obj.F(band);

X = gather(obj.GD3.X);
Y = gather(obj.GD3.Y);
Z = gather(obj.GD3.Z);


x_c = gather(x_shift);
y_c = gather(y_shift);
z_c = gather(z_shift);


tic;
%C(band) = interp3(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,C,x_shift,y_shift,z_shift,'linear');
C(band) = gpuArray(interp3(X,Y,Z,gather(C),x_c,y_c,z_c,'cubic'));
toc;

step = zeros(obj.GD3.Size, 'gpuArray');

mask = oldC < 0;
for i=1:100
	step = feval(obj.ENORK2_reinitiliaztion_step, step, C, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp = C - step;
%	Ctmp(band) = gpuArray(interp3(X,Y,Z,gather(Ctmp),x_c,y_c,z_c,'cubic'));
%	Ctmp(band) = interp3(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,Ctmp,x_shift,y_shift,z_shift,'linear');
	Ctmp = obj.ENORK2Extend(Ctmp,100);
	step = feval(obj.ENORK2_reinitiliaztion_step, step, Ctmp, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp2 = Ctmp - step;
%	Ctmp2(band) = gpuArray(interp3(X,Y,Z,gather(Ctmp2),x_c,y_c,z_c,'cubic'));
%	Ctmp2(band) = interp3(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,Ctmp2,x_shift,y_shift,z_shift,'linear');
	Ctmp2 = obj.ENORK2Extend(Ctmp2,100);
	C = (C + Ctmp2) / 2;

end

%C(band) = gpuArray(interp3(X,Y,Z,gather(C),x_c,y_c,z_c,'cubic'));
C = obj.ENORK2Extend(C, 100);

figure

subplot(2,2,1)
%obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6],obj.GD3.Z-0.3)
obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6],original)

subplot(2,2,2)
%obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6],oldC)
obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6],oldC)


subplot(2,2,3)
obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6],C)
%obj.plotIsoField([0,0.1,0.2,0.3,0.4,0.5,0.6],C)













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

























