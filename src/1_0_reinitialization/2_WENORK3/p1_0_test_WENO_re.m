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

figure(1)

tic
map.F = map.ENORK2Reinitialization(map.F, 100);
toc
map.GPUsetCalculusToolBox

subplot(2,2,1)
map.plotField(0,map.MeanCurvature-2/Radius)
subplot(2,2,2)
mask = abs(map.F) < 3*map.GD3.Ds;
map.MeanCurvature(~mask) = nan;
map.plotSurfaceField(map.MeanCurvature, 2/Radius, 1, 'red');



map.F = F;
tic
map.F = map.WENORK3Reinitialization(map.F, 100);
toc
map.GPUsetCalculusToolBox

subplot(2,2,3)
map.plotField(0,map.MeanCurvature-2/Radius)
subplot(2,2,4)
mask = abs(map.F) < 3*map.GD3.Ds;
map.MeanCurvature(~mask) = nan;
map.plotSurfaceField(map.MeanCurvature, 2/Radius, 1, 'red');

%fun = @(x,y,z) (sqrt(x.^2+z.^2+y.^2) - Radius);
%map.F = fun(x,y,z);
%map.setCalculusToolBox
%
%subplot(1,2,2)
%map.plotField(0,map.MeanCurvature)
%map.plotSurfaceField(map.MeanCurvature, 2/Radius, 1, 'red');





























