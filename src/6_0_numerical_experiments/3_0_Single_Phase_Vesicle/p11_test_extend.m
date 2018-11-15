% test extend scheme

r0 = .6;
xmin = -1.0; xmax = 1.0;

N = 64;
%N = 128;
xv = linspace(xmin,xmax,N);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);
x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);
grid = SD.GD3(x,y,z);

[azimuth,elevation,r] = cart2sph(x,y,z);

F = sqrt(x.^2 + y.^2 + z.^2) - r0;
map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.setCalculusToolBox4;

%C0 = sin(2.*azimuth).*sin(2.*elevation).*r;
%C_true = sin(2.*azimuth).*sin(2.*elevation).*r0;

C0 = z; % sin(elevation).*r
C_true = sin(elevation) .* r0;

%C_ex = map.WENO5RK3Extend(C0, 300);
C_ex = map.WENORK3Extend(C0, 300);
%G_lap = map.GD3.Laplacian4(G);
%G_Bilap = map.GD3.Laplacian4(G_lap);
%G_Bilap = map.WENO5RK3Extend(G_Bilap,100);

diff = C_ex - C_true;

E = sqrt(map.surfaceIntegral(diff.^2));
%E = sqrt(map.surfaceIntegral(G_Bilap.^2));
fprintf('N: %5d, E: %5.3e \n',N,E)

