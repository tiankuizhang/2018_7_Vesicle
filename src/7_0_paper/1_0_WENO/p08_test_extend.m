% test extend scheme
[EG1_16,EG2_16,EGM_16,EL1_16,EL2_16,ELM_16,EB1_16,EB2_16,EBM_16] = test_extend(16);
[EG1_32,EG2_32,EGM_32,EL1_32,EL2_32,ELM_32,EB1_32,EB2_32,EBM_32] = test_extend(32);
[EG1_64,EG2_64,EGM_64,EL1_64,EL2_64,ELM_64,EB1_64,EB2_64,EBM_64] = test_extend(64);
[EG1_28,EG2_28,EGM_28,EL1_28,EL2_28,ELM_28,EB1_28,EB2_28,EBM_28] = test_extend(128);

fprintf('\n convergence of the surface field \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EG1_16,EG2_16,EGM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EG1_32,log(EG1_16/EG1_32)/log(2),EG2_32,log(EG2_16/EG2_32)/log(2),EGM_32,log(EGM_16/EGM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EG1_64,log(EG1_32/EG1_64)/log(2),EG2_64,log(EG2_32/EG2_64)/log(2),EGM_64,log(EGM_32/EGM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EG1_28,log(EG1_64/EG1_28)/log(2),EG2_28,log(EG2_64/EG2_28)/log(2),EGM_28,log(EGM_64/EGM_28)/log(2))

fprintf('\n convergence of the surface field Laplacian \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EL1_16,EL2_16,ELM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EL1_32,log(EL1_16/EL1_32)/log(2),EL2_32,log(EL2_16/EL2_32)/log(2),ELM_32,log(ELM_16/ELM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EL1_64,log(EL1_32/EL1_64)/log(2),EL2_64,log(EL2_32/EL2_64)/log(2),ELM_64,log(ELM_32/ELM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EL1_28,log(EL1_64/EL1_28)/log(2),EL2_28,log(EL2_64/EL2_28)/log(2),ELM_28,log(ELM_64/ELM_28)/log(2))

fprintf('\n convergence of the surface field biLaplacian \n')
fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EB1_16,EB2_16,EBM_16)
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EB1_32,log(EB1_16/EB1_32)/log(2),EB2_32,log(EB2_16/EB2_32)/log(2),EBM_32,log(EBM_16/EBM_32)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EB1_64,log(EB1_32/EB1_64)/log(2),EB2_64,log(EB2_32/EB2_64)/log(2),EBM_64,log(EBM_32/EBM_64)/log(2))
fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EB1_28,log(EB1_64/EB1_28)/log(2),EB2_28,log(EB2_64/EB2_28)/log(2),EBM_28,log(EBM_64/EBM_28)/log(2))

% test extend scheme with a specific grid size
function [EG1,EG2,EGM,EL1,EL2,ELM,EB1,EB2,EBM] = test_extend(N)

	r0 = .5;
	xmin = -1.0; xmax = 1.0;
	
	xv = linspace(xmin,xmax,N);
	yv = xv;
	zv = xv;
	
	[x,y,z] = meshgrid(xv,yv,zv);
	x = gpuArray(x);
	y = gpuArray(y);
	z = gpuArray(z);
	grid = SD.GD3(x,y,z);
	
	% elevation ranges from pi/2 to -pi/2 (from north pole to south pole)
	% this is different from the inclination angle which ranges from 0 to pi
	[azimuth,elevation,r] = cart2sph(x,y,z);
	
	F = sqrt(x.^2 + y.^2 + z.^2) - r0;
	map = SD.SDF3(grid,x,y,z,F);
	map.setDistance
	map.setCalculusToolBox4;

	mask = abs(map.F) < 1.5*map.GD3.Ds;
	
	G_exact = sin(elevation) .* r0;
	G_lap_exact = -2/r0^2 .* G_exact;
	G_Bilap_exact = 4/r0^4 .* G_exact;
	
	G = map.WENO5RK3Extend(z,100);
	%G = map.WENORK3Extend(z,100);
	G_lap = map.GD3.Laplacian4(G);
	G_lap = map.WENO5RK3Extend(G_lap,100);
	G_Bilap = map.GD3.Laplacian4(G_lap);
	G_Bilap = map.WENO5RK3Extend(G_Bilap,100);
	
	% error in G
	DG = G - G_exact;
	DG = map.WENO5RK3Extend(DG,100);
	EG1 = map.surfaceIntegral(abs(DG));
	EG2 = sqrt(max(map.surfaceIntegral(DG.^2),0));
	EGM = max(abs(DG(mask)));
	%fprintf('%3d, EG1: %5.3e, EG2: %5.3e, EGM: %5.3e \n',N,EG1,EG2,EGM)

	% error in G_lap
	DL = G_lap - G_lap_exact;
	DL = map.WENO5RK3Extend(DL,100);
	EL1 = map.surfaceIntegral(abs(DL));
	EL2 = sqrt(max(map.surfaceIntegral(DL.^2),0));
	ELM = max(abs(DL(mask)));
%	fprintf('%3d, EL1: %5.3e, EL2: %5.3e, ELM: %5.3e \n',N,EL1,EL2,ELM)

	% error in G_Bilap
	DB = G_Bilap - G_Bilap_exact;
	DB = map.WENO5RK3Extend(DB,100);
	EB1 = map.surfaceIntegral(abs(DB));
	EB2 = sqrt(max(map.surfaceIntegral(DB.^2),0));
	EBM = max(abs(DB(mask)));
%	fprintf('%3d, EB1: %5.3e, EB2: %5.3e, EBM: %5.3e \n',N,EB1,EB2,EBM)

end
