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
Radius = 0.6;
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-Radius;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);


obj = map;
OldC = obj.GD3.Z;
NewC = OldC;

NewC = map.WENORK3Extend(OldC,100);

figure
iso = 0.3;
obj.plotSurface(0,0.8,'Green',1)
obj.plotSurfaceField(OldC,iso,0.8,'Red')
obj.plotSurfaceField(NewC,iso,0.8,'Blue')






