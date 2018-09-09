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
map.GPUsetCalculusToolBox;

A = z - 0.3;
Extend = map.ENORK2Extend(A,100);

tic
Re  = map.ENORK2Extend( map.ENORK2ClosetPointSurfaceRedistance(Extend,100,10), 100);
toc
tic
Sur = map.ENORK2Extend( map.ENORK2ClosetPointSurfaceRedistance(Extend,10,5), 100);
toc


map.A = Re;
map.AsetCalculusToolBox

subplot(2,2,1)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);
subplot(2,2,2)
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);

map.A = Sur;
map.AsetCalculusToolBox

subplot(2,2,3)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);
subplot(2,2,4)
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);




