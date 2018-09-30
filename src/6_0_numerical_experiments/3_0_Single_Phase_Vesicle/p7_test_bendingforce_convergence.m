% we shall test the convergence of various geometrical quantities
% calculated with FDM from a level set approach to see if there is
% any error in implementation

Size = [64,64,64];
%Size = [96,96,96];
%Size = [128,128,128];
%Size = [192,192,192];
%Size = [256,256,256];
[X,Y,Z,f,Nx,Ny,Nz,MC,SL_MC] = SD.Shape.RedBloodCell(Size,1.0);

grid = SD.GD3(X,Y,Z);
map = SD.SDF3(grid,X,Y,Z,f(X,Y,Z));
map.F = map.WENORK3Reinitialization(map.F,300);
map.setCalculusToolBox4

% theoretical Normal, MeanCurvature and Surface Laplacian of MeanCurvature
Nx = Nx(X,Y,Z);
Ny = Ny(X,Y,Z);
Nz = Nz(X,Y,Z);
MC = MC(X,Y,Z);
SL_MC = SL_MC(X,Y,Z);

mask = abs(map.F) < map.GD3.Ds * 2;
SL_MC = map.WENORK3Extend(SL_MC,100);


diffMC = map.WENORK3Extend(MC - map.MeanCurvature,100);

MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian4(map.WENORK3Extend(map.MeanCurvature,100));
diffSL_MC = map.WENORK3Extend(SL_MC - MeanCurvatureSurfaceLaplacian,100);

MCL2Error = sqrt(map.surfaceIntegral(diffMC.^2))
SLMCL2Error = sqrt(map.surfaceIntegral(diffSL_MC.^2))

