% test order of convergence of the new surface redistance scheme

r0 = 0.6;
xmin = -1.0; xmax = 1.0;
N = 128;
%N = 64;

xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);
[azimuth,elevation,r] = cart2sph(x,y,z);

F= sqrt(x.^2 + y.^2 + z.^2) - r0;
map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.setCalculusToolBox4;

map.A = z;
map.A = map.WENO5RK3Extend(map.A,200);
map.AsetCalculusToolBox4;

G_exact = elevation .* r0;
G_nu = map.WENO5RK3ClosetPointSurfaceRedistance(map.A, 300, 20);

diff = map.WENO5RK3Extend(G_exact-G_nu,100);

iso = -0.7:0.1:0.7;
%map.plotIsoField(iso,map.A,false)

E = sqrt(map.surfaceIntegral(diff.^2));
fprintf('N: %5d, E: %5.3e \n',N, E)

