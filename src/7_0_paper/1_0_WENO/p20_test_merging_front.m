[EF1_16,EF2_16,EFM_16,EMC1_16,EMC2_16,EMCM_16,EGC1_16,EGC2_16,EGCM_16,ELMC1_16,ELMC2_16,ELMCM_16] = CalculateError(16);
[EF1_32,EF2_32,EFM_32,EMC1_32,EMC2_32,EMCM_32,EGC1_32,EGC2_32,EGCM_32,ELMC1_32,ELMC2_32,ELMCM_32] = CalculateError(32);
[EF1_64,EF2_64,EFM_64,EMC1_64,EMC2_64,EMCM_64,EGC1_64,EGC2_64,EGCM_64,ELMC1_64,ELMC2_64,ELMCM_64] = CalculateError(64);
[EF1_28,EF2_28,EFM_28,EMC1_28,EMC2_28,EMCM_28,EGC1_28,EGC2_28,EGCM_28,ELMC1_28,ELMC2_28,ELMCM_28] = CalculateError(128);

fprintf('\n convergence of the level set function \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EF1_16,EF2_16,EFM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EF1_32,log(EF1_16/EF1_32)/log(2),EF2_32,log(EF2_16/EF2_32)/log(2),EFM_32,log(EFM_16/EFM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EF1_64,log(EF1_32/EF1_64)/log(2),EF2_64,log(EF2_32/EF2_64)/log(2),EFM_64,log(EFM_32/EFM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EF1_28,log(EF1_64/EF1_28)/log(2),EF2_28,log(EF2_64/EF2_28)/log(2),EFM_28,log(EFM_64/EFM_28)/log(2))

%fprintf('\n convergence of the Mean Curvature \n')
%fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
%fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EMC1_16,EMC2_16,EMCM_16)
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EMC1_32,log(EMC1_16/EMC1_32)/log(2),EMC2_32,log(EMC2_16/EMC2_32)/log(2),EMCM_32,log(EMCM_16/EMCM_32)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EMC1_64,log(EMC1_32/EMC1_64)/log(2),EMC2_64,log(EMC2_32/EMC2_64)/log(2),EMCM_64,log(EMCM_32/EMCM_64)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EMC1_28,log(EMC1_64/EMC1_28)/log(2),EMC2_28,log(EMC2_64/EMC2_28)/log(2),EMCM_28,log(EMCM_64/EMCM_28)/log(2))
%
%fprintf('\n convergence of the Gaussian Curvature \n')
%fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
%fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EGC1_16,EGC2_16,EGCM_16)
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EGC1_32,log(EGC1_16/EGC1_32)/log(2),EGC2_32,log(EGC2_16/EGC2_32)/log(2),EGCM_32,log(EGCM_16/EGCM_32)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EGC1_64,log(EGC1_32/EGC1_64)/log(2),EGC2_64,log(EGC2_32/EGC2_64)/log(2),EGCM_64,log(EGCM_32/EGCM_64)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EGC1_28,log(EGC1_64/EGC1_28)/log(2),EGC2_28,log(EGC2_64/EGC2_28)/log(2),EGCM_28,log(EGCM_64/EGCM_28)/log(2))
%
%fprintf('\n convergence of the surface Laplacian of the MeanCurvature \n')
%fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
%fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,ELMC1_16,ELMC2_16,ELMCM_16)
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,ELMC1_32,log(ELMC1_16/ELMC1_32)/log(2),ELMC2_32,log(ELMC2_16/ELMC2_32)/log(2),ELMCM_32,log(ELMCM_16/ELMCM_32)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,ELMC1_64,log(ELMC1_32/ELMC1_64)/log(2),ELMC2_64,log(ELMC2_32/ELMC2_64)/log(2),ELMCM_64,log(ELMCM_32/ELMCM_64)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,ELMC1_28,log(ELMC1_64/ELMC1_28)/log(2),ELMC2_28,log(ELMC2_64/ELMC2_28)/log(2),ELMCM_28,log(ELMCM_64/ELMCM_28)/log(2))
%
function [EF1,EF2,EFM,EMC1,EMC2,EMCM,EGC1,EGC2,EGCM,ELMC1,ELMC2,ELMCM] = CalculateError(N)
	Nx = N; Ny = N; Nz = 1.5*N;

	r = .6;
	xmin = -1.0; xmax = 1.0;
	xv = linspace(xmin,xmax,N);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;
	
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	grid = SD.GD3(x,y,z);
	
	z0 = 0.5 * r;
	F1 = sqrt( x.^2 + y.^2 + (z-z0).^2 ) - r; 
	F2 = sqrt( x.^2 + y.^2 + (z+z0).^2 ) - r; 
	F_exact = min(F1,F2);


	f1 = x.^2 + y.^2 + (z-z0).^2 - r^2;
	f2 = x.^2 + y.^2 + (z+z0).^2 - r^2;
	F = min(f1,f2);
	
	map = SD.SDF3(grid,x,y,z,F);
	map.setDistance
	map.F = map.WENO5RK3Reinitialization(map.F,1000);
	%map.F = map.WENORK3Reinitialization(map.F,1000);
	map.setCalculusToolBox4;

	mask = abs(map.F) < 1.5*map.GD3.Ds;
	
	% error in level set function
	DF = map.F - F_exact;
	DF = map.WENORK3Extend(DF,100);
	EF1 = map.surfaceIntegral(abs(DF));
	EF2 = sqrt(map.surfaceIntegral(DF.^2));
	EFM = max(abs(DF(mask)));
	
	EMC1=1; EMC2=1; EMCM=1;
	EGC1=1; EGC2=1; EGCM=1;
	ELMC1=1; ELMC2=1; ELMCM=1;
%	% error in Mean Curvature
%	MC = map.WENORK3Extend(map.MeanCurvature,100);
%	DMC = MC + 2.0 / r;
%	EMC1 = map.surfaceIntegral(abs(DMC));
%	EMC2 = sqrt(map.surfaceIntegral(DMC.^2));
%	EMCM = max(abs(DMC(mask)));
%
%	% error in Gaussian Curvature
%	GC = map.WENORK3Extend(map.GaussianCurvature,100);
%	DGC = GC - 1.0 / r.^2;
%	EGC1 = map.surfaceIntegral(abs(DGC));
%	EGC2 = sqrt(map.surfaceIntegral(DGC.^2));
%	EGCM = max(abs(DGC(mask)));
%
%	% error in SurfaceLaplacian of Mean curvature
%	SLMC = map.SurfaceLaplacian(MC);
%	%SLMC = map.SurfaceLaplacian4(MC);
%	SLMC = map.WENORK3Extend(SLMC,100);
%	ELMC1 = map.surfaceIntegral(abs(SLMC));
%	ELMC2 = sqrt(map.surfaceIntegral(SLMC.^2));
%	ELMCM = max(abs(SLMC(mask)));
%	
	%fprintf('%03d: MC Err: %5.3e, SLMC Err: %5.3e, F Err: %5.3e \n',N,E64,E64_2,EF)
	fprintf('%03d:    F1: %5.3e,    F2: %5.3e,    FM: %5.3e \n',N,EF1,EF2,EFM)
	%fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',N,EMC1,EMC2,EMCM)
	%fprintf('%03d:  EMC1: %5.3e,  EMC2: %5.3e,  EMCM: %5.3e \n',N,EGC1,EGC2,EGCM)
	%fprintf('%03d: ELMC1: %5.3e, ELMC2: %5.3e, ELMCM: %5.3e \n',N,ELMC1,ELMC2,ELMCM)

end
