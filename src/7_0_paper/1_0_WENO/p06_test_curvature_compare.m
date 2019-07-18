% draw iso surface of laplacian of mean curvature 
r = .6;
xmin = -1.0; xmax = 1.0;
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
%map.F = map.WENORK3Reinitialization(map.F,1000);
map.setCalculusToolBox4;
MC = map.WENORK3Extend(map.MeanCurvature,100);

SLMC = map.SurfaceLaplacian4(MC);
SLMC = map.WENORK3Extend(SLMC,100);

LMC = map.GD3.Laplacian4(map.MeanCurvature);
LMC_Ex = map.WENORK3Extend(LMC,100);

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
