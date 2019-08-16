
% we shall test the convergence of various geometrical quantities
% calculated with FDM from a level set approach to see if there is
% any error in implementation

Size16 = [16,16,16];
Size32 = [32,32,32];
Size64 = [64,64,64];
Size128 = [128,128,128];

%[~,~,~,f,Nx,Ny,Nz,MC,GC,SL_MC] = SD.Shape.RedBloodCell(Size16,1.0);

[EMC1_28,EMC2_28,EMCM_28,EGC1_28,EGC2_28,EGCM_28,ELMC1_28,ELMC2_28,ELMCM_28] = CalculateError(Size64,f,MC,GC,SL_MC);

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
	tmpSLMC = map.WENORK3Extend(SLMC,100);
	maxSLMC = max(abs(SLMC(mask)))
	MeanCurvatureSurfaceLaplacian = map.SurfaceLaplacian4(map.WENORK3Extend(map.MeanCurvature,100));
	DSLMC = map.WENORK3Extend(SLMC - MeanCurvatureSurfaceLaplacian,100);
	ELMC1 = map.surfaceIntegral(abs(DSLMC));
	ELMC2 = sqrt(map.surfaceIntegral(DSLMC.^2));
	ELMCM = max(abs(DSLMC(mask)));

%	fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',Size(1),EMC1,EMC2,EMCM)
%	fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',Size(1),EGC1,EGC2,EGCM)
%	fprintf('%03d: ELMC1: %5.3e, ELMC2: %5.3e, ELMCM: %5.3e \n',Size(1),ELMC1,ELMC2,ELMCM)
%	
	%map.plotField(0,DSLMC/maxSLMC,0.)
	map.plotSurface(0,1,'r','none')
	axis equal vis3d
	set(gca,'visible','off')
	%set(gca,'Color','k')
	%colorbar off
	title('red blood cell')
end













