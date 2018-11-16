% test extend scheme

r0 = .6;
xmin = -1.0; xmax = 1.0;

%N = 64;
%N = 128;
N = 196;
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

%C0 = sin(3.*azimuth).*sin(7.*elevation).*r;
%C_true = sin(3.*azimuth).*sin(7.*elevation).*r0;

%C0 = z; % sin(elevation).*r
%C_true = sin(elevation) .* r0;

%C_ex = map.WENO5RK3Extend(C0, 300);
%C_ex = map.WENORK3Extend(C0, 300);

G = map.WENO5RK3Extend(z,100);
%G = map.WENORK3Extend(z,100);
G_exact = sin(elevation) .* r0;

G_lap = map.GD3.Laplacian4(G);
G_lap = map.WENO5RK3Extend(G_lap,100);
G_lap_exact = -2/r0^2 .* G_exact;

G_Bilap = map.GD3.Laplacian4(G_lap);
G_Bilap = map.WENO5RK3Extend(G_Bilap,100);
G_Bilap_exact = 4/r0^4 .* G_exact;

%diff = C_ex - C_true;
%diff = G_lap - G_lap_exact;
diff = G_Bilap - G_Bilap_exact;


E = sqrt(map.surfaceIntegral(diff.^2));
%E = sqrt(map.surfaceIntegral(G_Bilap.^2));
fprintf('N: %5d, E: %5.3e \n',N,E)

