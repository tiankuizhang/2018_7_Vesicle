% verify order of accuracy for our scheme with Mean Curvature
r = .6;
xmin = -1.0; xmax = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 64;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,1000);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E64 = sqrt(map.surfaceIntegral(Diff.^2));
E64_2 = sqrt(map.surfaceIntegral(SLMC.^2));

DF = map.F - (sqrt(x.^2+y.^2+z.^2) - r);
DF = map.WENORK3Extend(DF,100);
mask = abs(map.F) < 2*map.GD3.Ds;
EF64 = sqrt(map.surfaceIntegral(DF.^2));

fprintf('64: MC Err: %5.3e, SLMC Err: %5.3e, F Err: %5.3e \n',E64,E64_2,EF64)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 96;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,1000);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E96 = sqrt(map.surfaceIntegral(Diff.^2));
E96_2 = sqrt(map.surfaceIntegral(SLMC.^2));

DF = map.F - (sqrt(x.^2+y.^2+z.^2) - r);
DF = map.WENORK3Extend(DF,100);
mask = abs(map.F) < 2*map.GD3.Ds;
EF96 = sqrt(map.surfaceIntegral(DF.^2));

fprintf('96: MC Err: %5.3e, SLMC Err: %5.3e, F Err: %5.3e \n',E96,E96_2,EF96)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 128;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,1000);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E128 = sqrt(map.surfaceIntegral(Diff.^2));
E128_2 = sqrt(map.surfaceIntegral(SLMC.^2));

DF = map.F - (sqrt(x.^2+y.^2+z.^2) - r);
DF = map.WENORK3Extend(DF,100);
mask = abs(map.F) < 2*map.GD3.Ds;
EF128 = sqrt(map.surfaceIntegral(DF.^2));

fprintf('128: MC Err: %5.3e, SLMC Err: %5.3e, F Err: %5.3e \n',E128,E128_2,EF128)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 192;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

F = x.^2 + y.^2 + z.^2 - r^2;

map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,1000);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

Diff = MC + 2.0 / r;
%figure(1);map.plotField(0,Diff)
E192 = sqrt(map.surfaceIntegral(Diff.^2));
E192_2 = sqrt(map.surfaceIntegral(SLMC.^2));

DF = map.F - (sqrt(x.^2+y.^2+z.^2) - r);
DF = map.WENORK3Extend(DF,100);
mask = abs(map.F) < 2*map.GD3.Ds;
EF192 = sqrt(map.surfaceIntegral(DF.^2));

fprintf('192: MC Err: %5.3e, SLMC Err: %5.3e, F Err: %5.3e \n',E192,E192_2,EF192)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n order of convergence of mean curvature L2 error ')
fprintf('\n\tE64/E96 is %5.3f, \n\tE96/E128 is %5.3f, \n\tE128/E192 is %5.3f', E64/E96, E96/E128, E128/E192)
fprintf('\n\tE64/E128 is %5.3f, \n\tE96/E192 is %5.3f \n',E64/E128,E96/E192)
