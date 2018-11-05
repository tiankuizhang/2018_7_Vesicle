% test new reinitialization scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 64;
xv = linspace(-1,1,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

r = 0.6;
F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.F = map.WENO5RK3Reinitialization(map.F,300);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E64 = sqrt(map.surfaceIntegral(Diff.^2))
E64_2 = sqrt(map.surfaceIntegral(SLMC.^2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 96;
xv = linspace(-1,1,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

r = 0.6;
F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.F = map.WENO5RK3Reinitialization(map.F,300);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E96 = sqrt(map.surfaceIntegral(Diff.^2))
E96_2 = sqrt(map.surfaceIntegral(SLMC.^2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 128;
xv = linspace(-1,1,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

r = 0.6;
F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.F = map.WENO5RK3Reinitialization(map.F,300);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E128 = sqrt(map.surfaceIntegral(Diff.^2))
E128_2 = sqrt(map.surfaceIntegral(SLMC.^2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E64/E96
E96/E128
E64/E128

E64_2
E96_2
E128_2
