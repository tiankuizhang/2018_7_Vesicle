% test ENO-RK2 reinitialization scheme

% create a 3D grid
xv = linspace(-1,1,256);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

% create a SDF3 instance
Radius = 0.6;
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-Radius;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%AF = map.GD3.Z-0.3;
%
%tic
%AF1 = map.ENORK2Extend(AF,0);
%toc
%
%tic
%AF2 = map.ENORK2CentralUpwindSurfaceRedistance(AF,100);
%toc
%
%tic
%AF3 = map.ENORK2CentralUpwindSurfaceRedistance(AF1,100);
%toc
%AF3 = map.ENORK2Extend(AF3,0);
%
%iso = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7];
%
%figure(1)
%subplot(1,2,1)
%map.plotIsoField(iso,AF)
%
%subplot(1,2,2)
%map.plotIsoField(iso,AF2)
%
%figure(2)
%subplot(1,2,1)
%map.plotIsoField(iso,AF1)
%
%subplot(1,2,2)
%map.plotIsoField(iso,AF3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check convergence of the surface redistance scheme

steady = @(azi,ele,ra) ele.*ra; % final steady state
[azimuth,elevation,r]=cart2sph(x,y,z); 

ThSurDisMap = steady(azimuth,elevation,r);

AF = map.ENORK2Extend(map.GD3.Z,200);
step = max(Radius*pi*0.5/(map.GD3.Dx*0.3), 100);
NuSurDisMap = map.ENORK2CentralUpwindSurfaceRedistance(AF,step);
NuSurDisMap = map.ENORK2Extend(NuSurDisMap,200);

iso = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7];

%figure(1)
%subplot(1,2,1)
%map.plotIsoField(iso,ThSurDisMap)
%subplot(1,2,2)
%map.plotIsoField(iso,NuSurDisMap)

mask = abs(map.F) < 2*map.GD3.Dx;
Ele = nnz(mask);

Diff = ThSurDisMap - NuSurDisMap;

L2Error = sqrt( sum( Diff(mask).^2 / Ele  )  );
MxError = max( map.GD3.Dx.^3.*abs( Diff(mask)  ) );

figure
map.plotField(0,Diff)












