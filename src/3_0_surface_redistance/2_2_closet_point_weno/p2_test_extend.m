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
A = z - 0.3;
map.A = A;

map.GPUsetCalculusToolBox;
map.AsetCalculusToolBox;

Geo1 = map.GeodesicCurvature;


figure

map.A = map.WENORK3Extend(A,100);
map.AsetCalculusToolBox;
Geo2 = map.GeodesicCurvature;

subplot(2,2,1)
map.plotField(0,Geo2-Geo1)
subplot(2,2,2)
map.plotField(0,map.A-A)

map.A = map.ENORK2Extend(A,100);
map.AsetCalculusToolBox;
Geo3 = map.GeodesicCurvature;

subplot(2,2,3)
map.plotField(0,Geo3-Geo1)
subplot(2,2,4)
map.plotField(0,map.A-A)
