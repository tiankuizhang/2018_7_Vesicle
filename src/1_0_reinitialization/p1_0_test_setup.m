% test setup of a SDF3 on a GPU

addpath(genpath('..'))

% create a 3D grid
xv = linspace(-5,5,64);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);



