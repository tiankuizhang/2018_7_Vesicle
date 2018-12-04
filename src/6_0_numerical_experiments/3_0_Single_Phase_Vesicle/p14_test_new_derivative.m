
r = .6;
xmin = -1.0; xmax = 1.0;

N = 64;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

%[x,y,z,f] = SD.Shape.Ellipsoid([128,128,64],0.65,"o");
grid = SD.GD3(x,y,z);
f = sqrt(x.^2+y.^2+z.^2) - r;
map = SD.SDF3(grid,x,y,z,f);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);

map.setCalculusToolBoxGA(0.04)
figure;map.plotField(0,map.GaussianCurvature-1/r^2)
%
%map.setCalculusToolBox
%figure;map.plotField(0,map.GaussianCurvature)
%
map.GPUsetCalculusToolBox
figure;map.plotField(0,map.GaussianCurvature-1/r^2)
%
map.setCalculusToolBox4
figure;map.plotField(0,map.GaussianCurvature-1/r^2)
