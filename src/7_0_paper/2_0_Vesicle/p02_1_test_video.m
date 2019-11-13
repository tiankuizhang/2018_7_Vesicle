%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu = SD.Simulation(mfilename, 'pear');
simu.simulationStart
pwd
Archive = true;
numFrame = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate single phase vesicle with constrained reduced volume and area difference

%TYPE = "o"; rd = 0.60; rad = .90; adrate = 0.00100; totalTime = 0.04; relaxTime = 0.04; % elliptocytes

%TYPE = "p"; rd = 0.55; rad = 1.6; adrate = 0.00100; totalTime = 1.4e-3; relaxTime = 1.05e-3; %necklaces 

%TYPE = "o"; rd = 0.55; rad = 1.3; adrate = 0.00100; totalTime = 0.002; relaxTime = 0.0004; % two leg star fish

%TYPE = "o"; rd = 0.55; rad = 1.3; adrate = 0.00100; totalTime = 0.04; relaxTime = 0.003; % four->two->three leg star fish

TYPE = "p"; rd = 0.90; rad = 1.1; adrate = 0.0001; totalTime = 0.016; relaxTime = 0.009; % pear

%TYPE = "o"; rd = 0.80; rad = .90; adrate = 0.00100; totalTime = 0.04; relaxTime = 0.032; % stomatocytes



GridSize = [64,64,64];
[x,y,z,f] = SD.Shape.Ellipsoid(GridSize,rd,TYPE,0.35);
grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,f);

map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.GPUsetCalculusToolBox

InitialArea = map.calArea;
EquivalentRadius = sqrt(InitialArea/(4*pi));

InitialVolume = map.calVolume;
ReducedVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);

MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
InitialAreaDifference = map.surfaceIntegral(MeanCurvature);
InitialReducedAreaDifference = - InitialAreaDifference/(8*pi*EquivalentRadius);

fprintf('area: %4.5f, volume: %4.5f, rv: %4.5f, rad: %4.5f\n', InitialArea, InitialVolume, ReducedVolume, InitialReducedAreaDifference)

FIG = figure('Name','Single Vesicle','Position',[10 10 1600 800])

textX = gather(map.GD3.xmin);
textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textZ = gather(map.GD3.zmin);

KappaB = 1.0; % bending rigidity
CFLNumber = .2;
filterWidth = gather(map.GD3.Ds)*5.0;

ExpectedAreaDifference = 8.*pi*EquivalentRadius * (-rad);

time = 0;
%frameTime = relaxTime;
frameTime = 0;
array_ene = [];
array_t = [];
i = 0;
%for i = 0:3000
while time < totalTime
	i = i+1;
	map.GPUsetCalculusToolBox

	z_shift = - (map.Box(5) + map.Box(6));
	y_shift = - (map.Box(3) + map.Box(4));
	x_shift = - (map.Box(1) + map.Box(2));

	% calculate geometry constraints
	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - InitialVolume) / InitialVolume;
	
	ReducedVolume = 100 * (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);

	map.MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	map.GaussianCurvature = map.WENORK3Extend(map.GaussianCurvature,100);
	MeanCurvature = map.MeanCurvature;
	GaussianCurvature = map.GaussianCurvature;
	CurrentAreaDifference = map.surfaceIntegral(MeanCurvature);
	CurrentReducedAreaDifference = - CurrentAreaDifference / (8*pi*EquivalentRadius); 

	% calculate bending forces
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(map.MeanCurvature); 
	MeanCurvatureSurfaceLaplacian = map.ENORK2Extend(MeanCurvatureSurfaceLaplacian,100);
	NormalSpeedBend = KappaB * (MeanCurvatureSurfaceLaplacian + ...
			0.5 * MeanCurvature.^3 - 2.0 * MeanCurvature .* GaussianCurvature);

	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;

	% now solve for Lagrange multipliers
	c11 = InitialArea;
	c12 = CurrentAreaDifference; c21 = c12;
	c13 = map.surfaceIntegral(GaussianCurvature); c31 = c13;
	c22 = map.surfaceIntegral(MeanCurvature.^2);
	c23 = map.surfaceIntegral(MeanCurvature.*GaussianCurvature); c32 = c23;
	c33 = map.surfaceIntegral(GaussianCurvature.^2);

	volumeChangeRate = (InitialVolume - CurrentVolume) / Dt;
	areaChangeRate = (InitialArea - CurrentArea) / Dt;

	%if i < 500
	if time<relaxTime
		areaDifferenceChangeRate = (InitialAreaDifference - CurrentAreaDifference) / (2*Dt);
	else
		tmp = ExpectedAreaDifference - CurrentAreaDifference;
		tmp = sign(tmp) * min(abs(tmp), adrate*abs(ExpectedAreaDifference) )/ 2.;
		areaDifferenceChangeRate = tmp / Dt;
	end

	s1 = volumeChangeRate + map.surfaceIntegral(NormalSpeedBend);
	s2 = - areaChangeRate + map.surfaceIntegral(NormalSpeedBend.*MeanCurvature);
	s3 = - areaDifferenceChangeRate + map.surfaceIntegral(NormalSpeedBend.*GaussianCurvature);

	PTA = [c11,c12,c13;c21,c22,c23;c31,c32,c33] \ [s1;s2;s3];
	Pressure = PTA(1); Tension = PTA(2); TensionDA = PTA(3);

	% now calculate normal speed
	normalSpeed = Tension .* MeanCurvature + Pressure + ...
				TensionDA .* GaussianCurvature - NormalSpeedBend;

	% time step the level set function
	%normalSpeedSmoothed = smoothGMRES(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = smoothFFT(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);
	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance
	
	ene = KappaB * c22;
	array_ene = [array_ene; ene];
	array_t = [array_t time];

	fprintf('iter: %5d, ene: %4.5f, ar: %+4.5f, vol: %+4.5f, rd: %4.5f, rad: %+4.5f\n', i, ene, DiffArea, DiffVolume, ReducedVolume, CurrentReducedAreaDifference)

	%if mod(i,20)==0 || i==2
	%if time>frameTime || i==2
	if time>frameTime
		%frameTime = frameTime + (totalTime - relaxTime)/numFrame;
		frameTime = frameTime + totalTime/numFrame;
		clf(FIG)

		subplot(2,2,[1,3])
		titlestr1 = [ sprintf('rd:%.2f,rad:%.2f,gr:(%d,%d,%d) \n P:%.2f,T:%.2f,TA:%.2f', ReducedVolume, CurrentReducedAreaDifference,GridSize(1),GridSize(2),GridSize(3),Pressure,Tension,TensionDA ) ];
		title(titlestr1)
		map.plotSurface(0,1,'green','none');
		%map.plotField(0,normalSpeedSmoothed,0.5)
		%map.plotField(0,map.AHeaviside,0.01)
		%map.plotField(0,Tension,0.01)
		%map.GD3.DrawBox

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
		titlestr2 = [ sprintf('shift:(%.2f,%.2f,%.2f)',y_shift,x_shift,z_shift ) ];
		title(titlestr2)

		subplot(2,2,4)
		area( array_t, array_ene )
		titlestr3 = [ sprintf('iter: %5d, time: %.3e, ene: %5.5f', i,time,ene ) ];
		title(titlestr3)

		drawnow

		if Archive 
			FIG.InvertHardcopy = 'off'; % preserve background color
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,10)==0
		map.F = circshift(map.F, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		map.setDistance
		map.F = map.WENO5RK3Reinitialization(map.F,200);
	end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu.simulationEnd
simu.processImage(10)
%SD.NE.processImage(10,'Elliptocyte')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt + Alpha * Dt * map.GD3.LBiLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

	[normalSpeedSmoothed,~,~,~,~] = gmres(Op, S, [], 1e-12, 300, @mfun);
	normalSpeedSmoothed = reshape(normalSpeedSmoothed, map.GD3.Size);

	% preconditioner
	function y = mfun(S)
		fftS = fftn(reshape(S,map.GD3.Size));
		fftS = fftS ./ (1 + Alpha * Dt * ...
				(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
		y = real(ifftn(fftS));
		y = reshape(y, [map.GD3.NumElt, 1]);
	end

end

% solve (Idt + alpha * Dt * BiLaplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothFFT(map, NormalSpeed, Dt, Alpha)

	fftS = fftn(NormalSpeed);
	fftS = fftS ./ (1 + Alpha * Dt * ...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2).^2 );
	normalSpeedSmoothed = real(ifftn(fftS));
	
end


% solve (Idt - alpha * Dt * Laplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothDiffusionGMRES(map, NormalSpeed, Dt, Alpha)

	% operator to be solved
	Op = map.GD3.Idt - Alpha * Dt * map.GD3.LLaplacian;
	% reshape RHS into a colum vector
	S = reshape(NormalSpeed, [map.GD3.NumElt, 1]);

	[normalSpeedSmoothed,~,~,~,~] = gmres(Op, S, [], 1e-12, 300, @mfun);
	normalSpeedSmoothed = reshape(normalSpeedSmoothed, map.GD3.Size);

	% preconditioner
	function y = mfun(S)
		fftS = fftn(reshape(S,map.GD3.Size));
		fftS = fftS ./ (1 + Alpha * Dt * ...
				(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2) );
		y = real(ifftn(fftS));
		y = reshape(y, [map.GD3.NumElt, 1]);
	end

end


% solve (Idt - alpha * Dt * Laplacian)^(-1) with GMRES preconditioned by FFT
function normalSpeedSmoothed = smoothDiffusionFFT(map, NormalSpeed, Dt, Alpha)

	fftS = fftn(NormalSpeed);
	fftS = fftS ./ (1 + Alpha * Dt * ...
			(map.GD3.kx.^2 + map.GD3.ky.^2 + map.GD3.kz.^2) );
	normalSpeedSmoothed = real(ifftn(fftS));
	
end



