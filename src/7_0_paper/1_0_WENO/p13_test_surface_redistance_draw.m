N = 64;

iter1 = 100;
iter2 = 20;

r0 = 0.6;
xmin = -1.0; xmax = 1.0;

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

%map.A = r0 .* sin(elevation);
A0 = 1 - exp(2*z);
map.A = A0;

G_exact = r0 .* elevation;
G_nu = map.WENO5RK3ClosetPointSurfaceRedistance(map.A, iter1, iter2);

iso = linspace(-10, 10, 100);

ax1 = subplot(1,2,1);
map.plotIsoField(iso, A0, true)
axis equal vis3d
set(gca,'Color','k')
%axisLim = [-0.10 0.10];
%set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
title('before')
colorbar off

ax2 = subplot(1,2,2);
map.plotIsoField(iso, G_nu, true)
axis equal vis3d
set(gca,'Color','k')
%axisLim = [-0.10 0.10];
%set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
title('after')
colorbar off
