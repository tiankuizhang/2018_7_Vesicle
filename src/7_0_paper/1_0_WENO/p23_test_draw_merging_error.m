r = 0.6;
NxList = [16,32,64,128];
%NxList = [16,32,64];

figure
for j = 1:length(NxList)
	[EF1,EF2,EFM] = CalculateError(j,NxList(j),r,0.5);
end



function [EF1,EF2,EFM] = CalculateError(j,N,r,ra)
	Nx = N; Ny = N; Nz = 2*N;

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
	
	z0 = ra * r;
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
	
	fprintf('%03d: F1: %5.3e, F2: %5.3e, FM: %5.3e, Di: %.3f \n',N,EF1,EF2,EFM,2*z0-2*r)

	clf
	map.plotField(0,log(abs(DF))/log(10),0)
	set(gca,'visible','off')
	drawnow
	saveas(gcf, [sprintf('%03d',Nx),'error'],'fig')

%	subplot(1,4,j)
%	map.plotField(0,log(abs(DF))/log(10),0)
%	%view(45,0)
%	set(gca,'visible','off')
%	drawnow
%	saveas(gcf, [sprintf('%03d',Nx),'error'],'fig')
end
