
[x,y,z,f] = SD.Shape.Ellipsoid([128,128,64],0.65,"o");
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,f);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);

