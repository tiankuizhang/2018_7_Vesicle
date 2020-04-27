% relaxation of ellispoid of different types and reduced volume

iteration = 1000;
Archive = false;
%GridSize = [48,48,48]; ReinitializationRate = 1;CFLNumber = 0.1;
GridSize = [64,64,64]; ReinitializationRate = 10;CFLNumber = 0.1;
%GridSize = [80,80,80]; ReinitializationRate = 1;CFLNumber = 0.1;
%GridSize = [96,96,96]; ReinitializationRate = 1;CFLNumber = 0.1;
%GridSize = [128,128,128];ReinitializationRate = 20;CFLNumber = 0.1;
%simu = none;

type = "P"; rv = 0.6; SponC = -20; SampleRate = 2;% iter: 400; totalTime: 1.8e-4

%totalTime = 2.5e-4;
%totalTime = 1.6e-3;
totalTime = 1.6e-2;
numFrame = 50;
SampleRate = 2;

% create the initial distance map
%[x,y,z,f] = SD.Shape.Ellipsoid([128,128,128],rv,type);
[x,y,z,f] = SD.Shape.Ellipsoid(GridSize,rv,type, 0.35);
load pear.mat
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,500);

map.GPUsetCalculusToolBox
InitialArea = map.calArea;
InitialVolume = map.calVolume;
ReduceVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);

fprintf('initial area: %4.5f, initial volume: %4.5f, reduced volume: %4.5f\n', ...
	InitialArea, InitialVolume, ReduceVolume)

% name and size of figure
FIG = figure('Name','Single Phase Vesicle','Position',[10 10 1600 800])

% position to show iteration number and accumulated time
textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

% bending mudulus
Kappa = 1.0;
filterWidth = gather(map.GD3.Ds)*5.0;

% dynamics
time = 0;
frameTime = totalTime/numFrame;
array_ene = [];
array_t = [];
array_da = [];
array_dv = [];
array_P = [];
array_ten = [];
i = 1;
while time < totalTime
%for i = 1:iteration
	i = i + 1;
	map.GPUsetCalculusToolBox

	z_shift = - (map.Box(5) + map.Box(6));
	y_shift = - (map.Box(3) + map.Box(4));
	x_shift = - (map.Box(1) + map.Box(2));

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - InitialVolume) / InitialVolume;
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);

	% extend mean curvature away from surface
	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);

	% surface Laplacian of mean curvature and numerical Hamiltonian
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(MeanCurvature); 
	NormalSpeedBend = Kappa * (MeanCurvatureSurfaceLaplacian ...
			+ 0.5 * MeanCurvature .* (MeanCurvature.^2 - 4 * map.GaussianCurvature) ...
			+ 2. * SponC * map.GaussianCurvature ...
			- 0.5 * SponC^2 * MeanCurvature );
	NormalSpeedBend = map.ENORK2Extend(NormalSpeedBend, 100);

	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	%Dt = 0.5e-5;
	time = time + Dt;


	% now solve for tension and pressure to enfore total area and volume
	c11 = map.surfaceIntegral(MeanCurvature.^2); % surface integral of mean curvature squared
	c12 = map.surfaceIntegral(MeanCurvature); % surface integral of mean curvature
	c21 = c12;
	c22 = InitialArea;

	%s1 = - map.surfaceIntegral(MeanCurvature .* NormalSpeedBend) ... % sign error
	s1 = map.surfaceIntegral(MeanCurvature .* NormalSpeedBend) ...
		 - (InitialArea - CurrentArea) / Dt;
	s2 = map.surfaceIntegral(NormalSpeedBend) + (InitialVolume - CurrentVolume) / Dt;
	
	TP = [c11,c12;c21,c22]\[s1;s2];
	Tension = TP(1);
	Pressure = TP(2);

	fprintf('iter: %5d, ene: %4.5f, ar: %4.5f, vol: %4.5f, rd: %4.5f\n', ...
			i, c11, DiffArea, DiffVolume, ReduceVolume)
	% now calculate normal Speed
	normalSpeed = (Tension * MeanCurvature - NormalSpeedBend + Pressure) .* map.FGradMag;
	normalSpeedSmoothed = map.GD3.smoothGMRES(normalSpeed, Dt, 0.5);
	%normalSpeedSmoothed = map.GD3.smoothFFT(normalSpeed, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);

	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance

	%ene = c11;
	ene = map.surfaceIntegral((MeanCurvature-SponC).^2); % surface integral of mean curvature squared
	if(i<=3), maxEne = gather(ene); end
	array_ene = [array_ene; ene];
	array_t = [array_t time];
	array_da = [array_da DiffArea/100];
	array_dv = [array_dv DiffVolume/100];
	array_P = [array_P Pressure];
	array_ten = [array_ten Tension];

	%if mod(i,SampleRate)==0 || time>totalTime
	%if mod(ceil(time/totalTime*numFrame),SampleRate)==0 || time>totalTime
	%if time>frameTime || i==2
	%if mod(i,SampleRate)==0 
	if mod(i,SampleRate)==0 
		frameTime = frameTime + totalTime/numFrame;
		timeStr = [sprintf('%04d: %0.5e, %0.5f', i,time,c11)];

		clf(FIG)

		subplot(2,2,[1,3])
		%titlestr1 = [ sprintf('rd:%.2f,rad:%.2f,gr:(%d,%d,%d) \n P:%.2f,T:%.2f,TA:%.2f', ReducedVolume, CurrentReducedAreaDifference,GridSize(1),GridSize(2),GridSize(3),Pressure,Tension,TensionDA ) ];
		titlestr1 = [ sprintf('rd:%.2f,gr:(%d,%d,%d) \n P:%.2f,T:%.2f', ReducedVolume, GridSize(1),GridSize(2),GridSize(3),Pressure,Tension ) ];
		title(titlestr1)
		map.plotSurface(0,1,'green','none');
		%map.plotField(0,normalSpeedSmoothed,0.5)
		%map.plotField(0,map.AHeaviside,0.01)
		%map.plotField(0,Tension,0.01)
		map.GD3.DrawBox

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		set(gca,'Color','k')

		zoom(1.0)

		subplot(2,2,2)
		xslice = ceil(map.GD3.ncols / 2);
		Fslice = reshape(map.F(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Y = reshape(map.GD3.Y(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Z = reshape(map.GD3.Z(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		contour(Y,Z,Fslice,[0,0],'blue','LineWidth',3)
		axis equal
		titlestr2 = [ sprintf('shift:(%.2f,%.2f,%.2f), SponC: %.2f',y_shift,x_shift,z_shift ,SponC) ];
		title(titlestr2)

		subplot(2,2,4)
		area( array_t, array_ene )
		ylim(gca,[0 maxEne])
		titlestr3 = [ sprintf('iter: %5d, time: %.3e, ene: %5.5f', i,time,ene ) ];
		title(titlestr3)


%		map.plotSurface(0,1,'Green','black')
%		%map.plotField(0,normalSpeedSmoothed,0.5)
%		ax = gca;
%		ax.Visible = 'off';
%		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
%		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')
%
%		zoom(1.0)
%
		drawnow

		if Archive 
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,ReinitializationRate)==0
		map.F = map.WENO5RK3Reinitialization(map.F,100);
	end

end





