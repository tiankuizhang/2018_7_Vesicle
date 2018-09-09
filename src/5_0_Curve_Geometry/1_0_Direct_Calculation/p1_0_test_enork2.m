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
Radius = 0.6;
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-Radius;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);
map.GPUsetCalculusToolBox;

A = z - 0.3;
Extend = map.ENORK2Extend(A,100);
Re = map.ENORK2Extend(map.ENORK2Reinitialization(Extend,100),100);
Sur = map.ENORK2Extend(map.ENORK2CentralUpwindSurfaceRedistance(Extend,100),100);

map.A = Extend;
map.AsetCalculusToolBox

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

subplot(2,2,1)
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);

subplot(2,2,2)
map.plotField(0,map.NormalCurvature)
map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);

subplot(2,2,3)
map.plotField(0,map.GeodesicTorsion)
map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);

subplot(2,2,4)
map.plotField(0,map.BPerpendicular)
map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare effects of extendsion,reinitialization and surface redistance scheme on
% calculation of curvature
%figure
%
%subplot(2,2,1)
%map.A = A;
%map.AsetCalculusToolBox
%map.plotField(0,map.GeodesicCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,2)
%map.A = Extend;
%map.AsetCalculusToolBox
%map.plotField(0,map.GeodesicCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,3)
%map.A = Re;
%map.AsetCalculusToolBox
%map.plotField(0,map.GeodesicCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,4)
%map.A = Sur;
%map.AsetCalculusToolBox
%map.plotField(0,map.GeodesicCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











