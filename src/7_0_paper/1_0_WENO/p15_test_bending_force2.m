
% we shall test the convergence of various geometrical quantities
% calculated with FDM from a level set approach to see if there is
% any error in implementation

Size16 = [16,16,16];
Size32 = [32,32,32];
Size64 = [64,64,64];
Size128 = [128,128,128];
%Size = [96,96,96];
%Size = [192,192,192];
%Size = [256,256,256];

%[~,~,~,f,Nx,Ny,Nz,MC,GC,SL_MC] = SD.Shape.RedBloodCell(Size16,1.0);

[EMC1_16,EMC2_16,EMCM_16,EGC1_16,EGC2_16,EGCM_16,ELMC1_16,ELMC2_16,ELMCM_16] = CalculateError(Size16, f,MC,GC,SL_MC);
[EMC1_32,EMC2_32,EMCM_32,EGC1_32,EGC2_32,EGCM_32,ELMC1_32,ELMC2_32,ELMCM_32] = CalculateError(Size32, f,MC,GC,SL_MC);
[EMC1_64,EMC2_64,EMCM_64,EGC1_64,EGC2_64,EGCM_64,ELMC1_64,ELMC2_64,ELMCM_64] = CalculateError(Size64, f,MC,GC,SL_MC);
[EMC1_28,EMC2_28,EMCM_28,EGC1_28,EGC2_28,EGCM_28,ELMC1_28,ELMC2_28,ELMCM_28] = CalculateError(Size128,f,MC,GC,SL_MC);

fprintf('\n convergence of the Mean Curvature \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EMC1_16,EMC2_16,EMCM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EMC1_32,log(EMC1_16/EMC1_32)/log(2),EMC2_32,log(EMC2_16/EMC2_32)/log(2),EMCM_32,log(EMCM_16/EMCM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EMC1_64,log(EMC1_32/EMC1_64)/log(2),EMC2_64,log(EMC2_32/EMC2_64)/log(2),EMCM_64,log(EMCM_32/EMCM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EMC1_28,log(EMC1_64/EMC1_28)/log(2),EMC2_28,log(EMC2_64/EMC2_28)/log(2),EMCM_28,log(EMCM_64/EMCM_28)/log(2))

fprintf('\n convergence of the Gaussian Curvature \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EGC1_16,EGC2_16,EGCM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EGC1_32,log(EGC1_16/EGC1_32)/log(2),EGC2_32,log(EGC2_16/EGC2_32)/log(2),EGCM_32,log(EGCM_16/EGCM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EGC1_64,log(EGC1_32/EGC1_64)/log(2),EGC2_64,log(EGC2_32/EGC2_64)/log(2),EGCM_64,log(EGCM_32/EGCM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EGC1_28,log(EGC1_64/EGC1_28)/log(2),EGC2_28,log(EGC2_64/EGC2_28)/log(2),EGCM_28,log(EGCM_64/EGCM_28)/log(2))

fprintf('\n convergence of the surface Laplacian of the MeanCurvature \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,ELMC1_16,ELMC2_16,ELMCM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,ELMC1_32,log(ELMC1_16/ELMC1_32)/log(2),ELMC2_32,log(ELMC2_16/ELMC2_32)/log(2),ELMCM_32,log(ELMCM_16/ELMCM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,ELMC1_64,log(ELMC1_32/ELMC1_64)/log(2),ELMC2_64,log(ELMC2_32/ELMC2_64)/log(2),ELMCM_64,log(ELMCM_32/ELMCM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,ELMC1_28,log(ELMC1_64/ELMC1_28)/log(2),ELMC2_28,log(ELMC2_64/ELMC2_28)/log(2),ELMCM_28,log(ELMCM_64/ELMCM_28)/log(2))


function [EMC1,EMC2,EMCM,EGC1,EGC2,EGCM,ELMC1,ELMC2,ELMCM] = CalculateError(Size,f,MC,GC,SL_MC)

	R = 1.0;

	% create a meshgrid
	Nx = Size(1);
	Ny = Size(2);
	Nz = Size(3);

	xmin = -1.3*R;
	xmax =  1.3*R;
	xv = linspace(xmin,xmax,Nx);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;

	[X,Y,Z] = meshgrid(xv,yv,zv);

	X = gpuArray(X);
	Y = gpuArray(Y);

	% initialize the distance map
	grid = SD.GD3(X,Y,Z);
	map = SD.SDF3(grid,X,Y,Z,f(X,Y,Z));
	map.setDistance

	map.F = map.WENORK3Reinitialization(map.F,1000);
	%map.F = map.WENO5RK3Reinitialization(map.F,1000);
	%map.F = map.WENORK3Reinitialization(map.F,300);
	map.setCalculusToolBox4

	mask = abs(map.F) < 1.5*map.GD3.Ds;

	% error in Mean Curvature
	MC = MC(X,Y,Z);
	DMC = map.WENORK3Extend(MC - map.MeanCurvature,100);
	EMC1 = map.surfaceIntegral(abs(DMC));
	EMC2 = sqrt(map.surfaceIntegral(DMC.^2));
	EMCM = max(abs(DMC(mask)));

	% error in Gaussian Curvature
	GC = GC(X,Y,Z);
	DGC = map.WENORK3Extend(GC - map.GaussianCurvature,100);
	EGC1 = map.surfaceIntegral(abs(DGC));
	EGC2 = sqrt(map.surfaceIntegral(DGC.^2));
	EGCM = max(abs(DGC(mask)));

	% error in SurfaceLaplacian of Mean curvature
	MC = map.WENORK3Extend(MC, 100);
	SLMC = map.SurfaceLaplacian4(MC);
	MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian4(map.WENORK3Extend(map.MeanCurvature,100));
	DSLMC = map.WENORK3Extend(SLMC - MeanCurvatureSurfaceLaplacian,100);
	ELMC1 = map.surfaceIntegral(abs(DSLMC));
	ELMC2 = sqrt(map.surfaceIntegral(DSLMC.^2));
	ELMCM = max(abs(DSLMC(mask)));

%	fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',Size(1),EMC1,EMC2,EMCM)
%	fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',Size(1),EGC1,EGC2,EGCM)
%	fprintf('%03d: ELMC1: %5.3e, ELMC2: %5.3e, ELMCM: %5.3e \n',Size(1),ELMC1,ELMC2,ELMCM)
%	
end













