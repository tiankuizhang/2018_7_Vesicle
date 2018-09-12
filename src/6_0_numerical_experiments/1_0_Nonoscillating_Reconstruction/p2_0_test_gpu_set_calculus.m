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


