% test extend scheme from A field 

% create the initial distance map
[x,y,z,f] = SD.Shape.Ellipsoid([64,64,128],0.85,"p");
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,f);
map.A = z;
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,100,50);

map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox

iso = gather(-20*map.GD3.Ds:4*map.GD3.Ds:20*map.GD3.Ds);
%map.plotIsoField(iso,map.A,true)
%map.plotField(0,map.GD3.X,0.01)

[az,ele,r] = cart2sph(map.GD3.X,map.GD3.Y,map.GD3.Z);

C = map.GD3.Z + map.GD3.X.^2 + map.GD3.Y.^3;
NewC = map.AENORK2Extend(C,50,200,50);

figure
subplot(1,3,1)
%map.plotField(0,sign(map.A).*map.nx,0.01)
map.plotField(0,C,0.01)
subplot(1,3,2)
%map.plotField(0,sign(map.A).*map.ny,0.01)
map.plotField(0,NewC,0.01)
subplot(1,3,3)
map.plotField(0,C-NewC,0.01)






