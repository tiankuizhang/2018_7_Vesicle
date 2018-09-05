% test ENO-RK2 reinitialization scheme

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
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-0.6;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);


%%%%

AF = map.GD3.Z-0.3;

tic
NewAF = map.ENORK2CentralUpwindSurfaceRedistance(AF,100);
toc

iso = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7];

figure
subplot(1,2,1)
map.plotIsoField(iso,AF)

subplot(1,2,2)
map.plotIsoField(iso,NewAF)















