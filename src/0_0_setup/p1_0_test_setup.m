% test setup of a SDF3 on a GPU

% addpath(genpath('..'))

% create a 3D grid
xv = linspace(-5,5,64);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

% create a SDF3 instance
value = sqrt(x.^2+y.^2+z.^2) - 3;
map = SD.SDF3(grid, x, y, z, value);




