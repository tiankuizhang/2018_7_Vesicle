
	N = 128;

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
	
	G = map.WENO5RK3Extend(z,500);
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

	% plot before and after
	FIG = figure('Name','Extend','Position',[10 10 1600 800]);

	ax1 = subplot(1,2,1);
	hold on
	map.plotField(0,z,0)
	map.plotSurfaceField(z,0,0.5,'r')
	map.plotSurfaceField(z,0.2,0.5,'r')
	map.plotSurfaceField(z,0.4,0.5,'r')
	map.plotSurfaceField(z,-0.2,0.5,'r')
	map.plotSurfaceField(z,-0.4,0.5,'r')
	axis equal vis3d
	set(gca,'Color','k')
	colorbar off
	title('before')
	axisLim = [-1.0 1.0];
	set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
	hold off
	
	ax2 = subplot(1,2,2);
	hold on
	map.plotField(0,G,0)
	map.plotSurfaceField(G,0,0.5,'r')
	map.plotSurfaceField(G,0.2,0.5,'r')
	map.plotSurfaceField(G,0.4,0.5,'r')
	map.plotSurfaceField(G,-0.2,0.5,'r')
	map.plotSurfaceField(G,-0.4,0.5,'r')
	axis equal vis3d
	set(gca,'Color','k')
	colorbar off
	title('after')
	axisLim = [-1.0 1.0];
	set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
	hold off

%	set(ax1, 'Position', [-0.2 0 0.7 0.7])
%	set(ax2, 'Position', [0.2 0 0.7 0.7])




