% test order of convergence of the new surface redistance scheme
%[EG1_16,EG2_16,EGM_16,Ec1_16,Ec2_16,EcM_16] = test_surface_redistance(16);
%[EG1_32,EG2_32,EGM_32,Ec1_32,Ec2_32,EcM_32] = test_surface_redistance(32);
%[EG1_64,EG2_64,EGM_64,Ec1_64,Ec2_64,EcM_64] = test_surface_redistance(64);
[EG1_28,EG2_28,EGM_28,Ec1_28,Ec2_28,EcM_28] = test_surface_redistance(128);
%
%fprintf('\n convergence of the geodesics \n')
%fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
%fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,EG1_16,EG2_16,EGM_16)
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,EG1_32,log(EG1_16/EG1_32)/log(2),EG2_32,log(EG2_16/EG2_32)/log(2),EGM_32,log(EGM_16/EGM_32)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,EG1_64,log(EG1_32/EG1_64)/log(2),EG2_64,log(EG2_32/EG2_64)/log(2),EGM_64,log(EGM_32/EGM_64)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,EG1_28,log(EG1_64/EG1_28)/log(2),EG2_28,log(EG2_64/EG2_28)/log(2),EGM_28,log(EGM_64/EGM_28)/log(2))
%
%fprintf('\n convergence of the geodesic curvature \n')
%fprintf('\t \t L1 error \t order \t L2 error \t order \t LM error \t order \n')
%fprintf('\t %03d \t %5.3e \t - \t %5.3e \t - \t %5.3e \t - \n',16,Ec1_16,Ec2_16,EcM_16)
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 32,Ec1_32,log(Ec1_16/Ec1_32)/log(2),Ec2_32,log(Ec2_16/Ec2_32)/log(2),EcM_32,log(EcM_16/EcM_32)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n', 64,Ec1_64,log(Ec1_32/Ec1_64)/log(2),Ec2_64,log(Ec2_32/Ec2_64)/log(2),EcM_64,log(EcM_32/EcM_64)/log(2))
%fprintf('\t %03d \t %5.3e \t %3.2f \t %5.3e \t %3.2f \t %5.3e \t %3.2f \n',128,Ec1_28,log(Ec1_64/Ec1_28)/log(2),Ec2_28,log(Ec2_64/Ec2_28)/log(2),EcM_28,log(EcM_64/EcM_28)/log(2))



function [EG1,EG2,EGM,Ec1,Ec2,EcM] = test_surface_redistance(N)
	iter1 = 400;
	iter2 = 100;

	r0 = 0.6;
	xmin = -1.0; xmax = 1.0;
	ele = pi/12;
	alpha = pi/4;
	%ele = 0;
	z0 = r0 * sin(ele);

	
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
	%map.A = exp(z) - 1;
	%map.A = exp(2*z) - exp(2* r0*sin(ele) );
	%map.A = exp(2*(z-z0)) - 1.0;
	%map.A = z-z0;
	%map.A = exp( 2.0*r0*(elevation - ele) ) - 1;
	%map.A = exp(2*z) - 1;
	map.A = elevation - ele;
	
	G_exact = r0 .* (elevation - ele);
	G_nu = map.WENO5RK3ClosetPointSurfaceRedistance(map.A, iter1, iter2);
	%G_nu = map.WENORK3ClosetPointSurfaceRedistance(map.A, 100, 100);
	%G_nu = map.ENORK2ClosetPointSurfaceRedistance(map.A, 100, 100);
	map.A = G_nu;
	map.AsetCalculusToolBox4;
	
	GeodesicCurvature_exact = tan(elevation) ./ r0;

	%mask = abs(GeodesicCurvature_exact) > 2;
	mask = abs(elevation-ele) > alpha; % region to be excluded
	mask2 = abs(map.F) < 1.5*map.GD3.Ds;
	mask3 = abs(map.A) < 1.5*map.GD3.Ds;

	% error in G
	DG = G_nu - G_exact;
	DG = map.WENO5RK3Extend(DG,100);
%	DG = map.WENORK3Extend(DG,100);
	DG(mask) = 0.0;
	EG1 = map.surfaceIntegral(abs(DG));
	EG2 = sqrt(map.surfaceIntegral(DG.^2));
	EGM = max(abs(DG(~mask&mask2)));

%	EG1 = map.LineIntegral(abs(DG));
%	EG2 = sqrt(map.LineIntegral(DG.^2));
%	EGM = max(abs(DG(mask2&mask3)));
%
	% error in curvature
	Dc = GeodesicCurvature_exact - map.GeodesicCurvature;
	Dc = map.WENO5RK3Extend(Dc,100);
	Dc(mask) = 0.0;
	Ec1 = map.surfaceIntegral(abs(Dc));
	Ec2 = sqrt(map.surfaceIntegral(Dc.^2));
	EcM = max(abs(Dc(~mask&mask2)));

%	Ec1 = map.LineIntegral(abs(Dc));
%	Ec2 = sqrt(map.LineIntegral(Dc.^2));
%	EcM = max(abs(Dc(mask2&mask3)));

	%fprintf('N: %5d, EG1: %5.3e, EG2: %5.3e, EGM: %5.3e \n',N, EG1, EGl,EGM)
	fprintf('N: %5d, EG1: %5.3e, EG2: %5.3e, EGM: %5.3e, ',N, EG1, EG2,EGM)
	fprintf('Ec1: %5.3e, Ec2: %5.3e, EcM: %5.3e \n', Ec1, Ec2, EcM)

	if(N==128)
		fprintf('iter1: %3d, iter2: %3d, ele: %d, alpha: %d \n',iter1,iter2,pi/ele,pi/alpha)
	end
end

