% create a 3D grid
xv = linspace(-5,5,128);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

% create a SDF3 instance
Radius = 3.;

%fun = @(x,y,z) (0.1+(x-3.5).^2+(sqrt(y.^2+z.^2)-2).^2) .* (sqrt(x.^2/4+(z.^2+y.^2)/9)-1);
%fun = @(x,y,z) (0.1+(x-3.5).^2+(sqrt(y.^2+z.^2)-2).^2) .* (sqrt(x.^2+z.^2+y.^2) - Radius);
fun = @(x,y,z) x.^2+z.^2+y.^2 - Radius^2;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);
map.A = z;

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox

%figure(1)
%
%subplot(1,2,1)
%map.plot
%
%map.F = map.WENORK3Reinitialization(map.F, 100);
%map.GPUsetCalculusToolBox
%map.GPUAsetCalculusToolBox
%
%subplot(1,2,2)
%map.plot
%
%figure(2)
%subplot(1,2,1)
%map.plotField(0,map.MeanCurvature)
%map.plotSurfaceField(map.MeanCurvature, 2/Radius, 1, 'red');
%
%fun = @(x,y,z) (sqrt(x.^2+z.^2+y.^2) - Radius);
%map.F = fun(x,y,z);
%
%map.setCalculusToolBox
%subplot(1,2,2)
%map.plotField(0,map.MeanCurvature)
%%map.plotSurfaceField(map.MeanCurvature, 2/Radius, 1, 'red');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	fun = @(x,y,z) (sqrt(x.^2+z.^2+y.^2) - Radius);
	map.F = fun(x,y,z);
	map.F = x;

	obj = map;
	Fgpu = map.F;

	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

	v0 = zeros(obj.GD3.Size, 'gpuArray');
	v1 = zeros(obj.GD3.Size, 'gpuArray');
	v2 = zeros(obj.GD3.Size, 'gpuArray');
	v3 = zeros(obj.GD3.Size, 'gpuArray');

	[v0, v1, v2, v3, xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.cubic_boundary_correction, ...
			v0, v1, v2, v3, ...
			xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	


	nxpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	nxpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	nypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	nypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	nzpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	nzpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

	[nxpr, nxpl, nypf, nypb, nzpu, nzpd] = feval(obj.ENORK2_boundary_correction, ...
			nxpr, nxpl, nypf, nypb, nzpu, nzpd, Fgpu, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

	mask = nxpr ~= obj.GD3.Dx;
	boundary = obj.GD3.ooo(mask);
	current = boundary(1)

	[v0(current), v1(current), v2(current), v3(current); ...
	 ypf(current), ypb(current),zpu(current),zpd(current)]

	 [xpr(current), xpl(current); ...
	 nxpr(current),nxpl(current)]
	

	diff = xpr - nxpr;
	max(abs(diff(:)))



























