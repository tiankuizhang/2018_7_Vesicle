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
Extend = map.WENORK3Extend(A,100);

tic
tmp  = map.WENORK3ClosetPointSurfaceRedistance(Extend,100,10);
toc
Re  = map.WENORK3Extend(tmp, 100);

tic
tmp  = map.WENORK3ClosetPointSurfaceRedistance(Extend,5,5);
toc
Sur  = map.WENORK3Extend(tmp, 100);


iso = -0.7:0.1:0.7;
map.A = Re;
map.AsetCalculusToolBox

figure(1)
subplot(2,2,1)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);

subplot(2,2,2)
map.GeodesicCurvature = sign(map.GeodesicCurvature) .* min(abs(map.GeodesicCurvature),5);
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);

map.A = Sur;
map.AsetCalculusToolBox

subplot(2,2,3)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);

subplot(2,2,4)
map.GeodesicCurvature = sign(map.GeodesicCurvature) .* min(abs(map.GeodesicCurvature),5);
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);




