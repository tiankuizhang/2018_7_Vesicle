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

N = 100;

tic
for i=1:N
	map.setCalculusToolBox;
end
toc

tic
for i=1:N
	map.AsetCalculusToolBox;
end
toc

obj = SD.SDF3(grid, x, y, z, F);
obj.A = z;

tic
for i=1:N
	obj.GPUsetCalculusToolBox;
end
toc
tic
for i=1:N
	obj.GPUAsetCalculusToolBox;
end
toc

