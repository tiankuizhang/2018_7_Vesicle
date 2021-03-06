
% test new scheme to account for protein dependent properties for single phase vesicle
% reduced volume is fixed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load pear.mat
%simu = SD.Simulation(mfilename, 'pear_protein_pinch');
%simu.simulationStart
Archived = false
%pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters
iteration = 2000; relaxIter = 0;
GridSize = [64,64,64]; 
Kappa0 = 1; Kappa1 = 0.0; % bending modulus
%C0 = -0.0; C1 = -1.0; proteinCoverage = 1.0;
C0 = -10; C1 = .0; proteinCoverage = .0;
%Mu = 1000; % incompressibility of vesicle
Mu = 1000; % incompressibility of vesicle
CFLNumber = 0.2;
MinimumTimeStep = 1 * 1e-6; % for downward pear
%MinimumTimeStep = 1 * 1e-5; % for upward pear
MinimumTimeStep2 = 1 * 1e-6; % minimum time step after adding protein
%MinimumTimeStep2 = 0.00;
RelativeTimeScale = .1; % relative drag coefficient for protein motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ReducedVolume0 = 0.90; VesicleTYPE = "p"; ratio = 0.35;
%ReducedVolume0 = 0.90; VesicleTYPE = "p"; ratio = 0.35; PresetP = 162.51; PresetTDA = 19.71; ConserveRAD = false;
%ReducedVolume0 = 0.75; VesicleTYPE = "p"; ratio = 0.35; PresetP = 700; ConserveRAD = true;
%ReducedVolume0 = 0.8; VesicleTYPE = "o"; ratio = 0.2; ConserveRAD = false; C0 = 13;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
[x,y,z,~] = SD.Shape.Ellipsoid(GridSize,ReducedVolume0,VesicleTYPE,ratio);
Grid = SD.GD3(x,y,z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map = SD.SDF3(Grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.GPUsetCalculusToolBox

InitialArea = map.calArea;
EquivalentRadius = sqrt(InitialArea/(4*pi));
InitialVolume = map.calVolume;
InitialReducedVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);
expectedVolume = InitialVolume; % no water crossing membrane

MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
InitialAreaDifference = map.surfaceIntegral(MeanCurvature);
InitialReducedAreaDifference = - InitialAreaDifference/(8*pi*EquivalentRadius);

fprintf('initial area: %.3e, initial volume: %.3e, rd: %.3e\n', InitialArea, InitialVolume, InitialReducedVolume)
FIG = figure('Name','Vesicle','Position',[10 10 1600 800])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration
time = 0;
array_t = []; % time array
array_b = []; % bending energy
array_c = []; % area compressibility

localArea = ones(map.GD3.Size,'gpuArray');
protein = ones(map.GD3.Size,'gpuArray');
%protein = proteinCoverage * ones(map.GD3.Size,'gpuArray');
for i = 0:iteration
	map.GPUsetCalculusToolBox

	z_shift = - (map.Box(5) + map.Box(6));
	y_shift = - (map.Box(3) + map.Box(4));
	x_shift = - (map.Box(1) + map.Box(2));

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - expectedVolume) / expectedVolume;
	ReducedVolume = (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);

	Kappa = Kappa0 + Kappa1 .* protein; % bending modulus field
	SC = C0 + C1 .* protein; % spontaneous curvature field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now calculate normal and tangential forces
	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	GaussianCurvature = map.WENORK3Extend(map.GaussianCurvature,100);
	localTension = 0.5*Mu .* ... % due to incompressibility
		( (1.0./localArea - 1).*(3.0./localArea - 1) + (1 - localArea.^2) );
	proTension = - 0.5 * Kappa1 * protein .* (MeanCurvature - SC).^2 ...
					+ C1 * protein .* Kappa .* (MeanCurvature - SC);

	CurrentAreaDifference = map.surfaceIntegral(MeanCurvature);
	CurrentReducedAreaDifference = - CurrentAreaDifference / (8*pi*EquivalentRadius); 
	% (minus) bending forces
	NormalBendSpeed = map.GD3.Laplacian( Kappa.*(MeanCurvature-SC) ) + ...
		+ Kappa .* (0.5*MeanCurvature.^3 - 2*MeanCurvature.*GaussianCurvature ...
					 + 2*SC.*GaussianCurvature - 0.5*SC.^2.*MeanCurvature) ...
		- MeanCurvature .* (localTension + proTension) ;
	NormalBendSpeed = map.ENORK2Extend(NormalBendSpeed,100);
	% tangential forces
	[tempx1,tempy1,tempz1] = map.GD3.Gradient(Kappa);
	[tempx2,tempy2,tempz2] = map.GD3.Gradient(Kappa.*SC);
	[tempx3,tempy3,tempz3] = map.GD3.Gradient( ...
			0.5*Kappa.*SC.^2 + localTension + proTension);
	tvx = 0.5 * MeanCurvature.^2 .* tempx1 - MeanCurvature .* tempx2 + tempx3;
	tvy = 0.5 * MeanCurvature.^2 .* tempy1 - MeanCurvature .* tempy2 + tempy3;
	tvz = 0.5 * MeanCurvature.^2 .* tempz1 - MeanCurvature .* tempz2 + tempz3;
	tv = sqrt(tvx.^2 + tvy.^2 + tvz.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select time step
	mask = abs(map.F) < 2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalBendSpeed(mask)));
	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	Dt = max(MinimumTimeStep,Dt)
	time = time + Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contrain total area and volume
%	c11 = map.surfaceIntegral(MeanCurvature.^2);
%	c12 = map.surfaceIntegral(MeanCurvature); c21 = c12;
%	c22 = CurrentArea;
%
%	areaChangeRate = (InitialArea - CurrentArea) / Dt;
%	volumeChangeRate = (InitialVolume - CurrentVolume) / Dt;
%	s1 = - areaChangeRate + map.surfaceIntegral(NormalBendSpeed.*MeanCurvature);
%	s2 = volumeChangeRate + map.surfaceIntegral(NormalBendSpeed);
%
%	TP = [c11,c12;c21,c22] \ [s1;s2];
%	Tension = TP(1); Pressure= TP(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% more complex constraints
	c11 = CurrentArea;
	c12 = CurrentAreaDifference; c21 = c12;
	c13 = map.surfaceIntegral(GaussianCurvature); c31 = c13;
	c22 = map.surfaceIntegral(MeanCurvature.^2);
	c23 = map.surfaceIntegral(MeanCurvature.*GaussianCurvature); c32 = c23;
	c33 = map.surfaceIntegral(GaussianCurvature.^2);

	volumeChangeRate = (InitialVolume - CurrentVolume) / Dt;
	areaChangeRate = (InitialArea - CurrentArea) / Dt;
	areaDifferenceChangeRate = (InitialAreaDifference - CurrentAreaDifference) / (2*Dt);

	s1 = volumeChangeRate + map.surfaceIntegral(NormalBendSpeed);
	s2 = - areaChangeRate + map.surfaceIntegral(NormalBendSpeed.*MeanCurvature);
	s3 = - areaDifferenceChangeRate + map.surfaceIntegral(NormalBendSpeed.*GaussianCurvature);

	% conserve volume, area and area difference
	%PTA = [c11,c12,c13;c21,c22,c23;c31,c32,c33] \ [s1;s2;s3];
	%Pressure = PTA(1); Tension = PTA(2); TensionDA = PTA(3);

%	if ConserveRAD
%		% conserve area and area difference, set volume to be constant
%		Pressure = PresetP;
%		TA = [c22,c23;c32,c33] \ [s2-Pressure*c21; s3-Pressure*c31];
%		Tension = TA(1); TensionDA = TA(2);
%	else
%		% conserve only area, set pressure to be constant
%		Pressure = PresetP; TensionDA =  PresetTDA;
%		Tension = (s2 - Pressure*c21 - TensionDA*c23) / c22; 
%	end

	PT = [c11,c12;c21,c22] \ [s1;s2];
	Pressure = PT(1); Tension = PT(2); TensionDA = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the complete normal speed 
	normalSpeed = Tension .* MeanCurvature + Pressure  ...
					+ TensionDA .* GaussianCurvature - NormalBendSpeed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for the level set function
	%maxKappa = max(abs(Kappa(mask)))
	maxKappa = 1; 
	%levelSetTimeStep = normalSpeed.*map.FGradMag - 1.0 * map.FRegularization(true);
	levelSetTimeStep = normalSpeed.*map.FGradMag ;
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, maxKappa);
	% do not forget to extend this term
	levelSetTimeStep = map.ENORK2Extend(levelSetTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divergence of the flow field
	%normalSpeedSmoothed = map.GD3.smoothFFT(normalSpeed, Dt, maxKappa);
	%Divergence =  - MeanCurvature .* normalSpeedSmoothed ...
	Divergence =  - MeanCurvature .* normalSpeed ...
				+ map.GD3.GradDotGrad(0.5*MeanCurvature.^2, Kappa) ...
				+ 0.5*MeanCurvature.^2 .* map.GD3.Laplacian(Kappa) ...
				- map.GD3.GradDotGrad(MeanCurvature, Kappa.*SC) ...
				- MeanCurvature .* map.GD3.Laplacian(Kappa.*SC) ...
				+ map.GD3.Laplacian(0.5*Kappa.*SC.^2 + localTension + proTension);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for localArea, i.e., sum of numerical fluxs
	laFlux = localArea .* Divergence + map.GD3.WENODotGrad(tvx,tvy,tvz,localArea);
	%maxTv = max(abs(tv(mask)));
	maxTv = 1.0;
	localAreaTimeStep = map.GD3.smoothDiffusionFFT(laFlux, Dt, maxTv);
	localAreaTimeStep = map.WENORK3Extend(localAreaTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divergence of protein motion field
	proDivergence = Divergence +  RelativeTimeScale * map.GD3.Laplacian(proTension) ;
% additional tangential motion
	[ptvx,ptvy,ptvz] = map.GD3.Gradient(RelativeTimeScale * proTension);
%	proDivergence = Divergence; ptvx = 0; ptvy = 0; ptvz = 0;
% (minus) time rate of change for protein field
	proFlux = protein .* proDivergence ...
		+ map.GD3.WENODotGrad(tvx+ptvx, tvy+ptvy, tvz+ptvz, protein);
	pTv = sqrt( (tvx+ptvx).^2 + (tvy+ptvy).^2 + (tvz+ptvz).^2 );
	%pmaxTv = max(abs(pTv(mask)));
	pmaxTv = 1.;
	proteinTimeStep = map.GD3.smoothDiffusionFFT(proFlux, Dt, pmaxTv);
	proteinTimeStep = map.WENORK3Extend(proteinTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divergence of protein motion field, only tangential effects are accounted for
%	proDivergence = map.GD3.Laplacian(proTension);
% additional tangential motion
%	[ptvx,ptvy,ptvz] = map.GD3.Gradient(proTension);
% (minus) time rate of change for protein field
%	proFlux = protein .* proDivergence ...
%		+ map.GD3.DotGrad(ptvx, ptvy, ptvz, protein);
%	proteinTimeStep = map.GD3.smoothDiffusionFFT(proFlux, Dt, maxTv);
%	proteinTimeStep = map.WENORK3Extend(proteinTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time step the system
	if i == relaxIter
		localArea = ones(map.GD3.Size,'gpuArray');
		protein = proteinCoverage * ones(map.GD3.Size,'gpuArray');
		MinimumTimeStep = MinimumTimeStep2;
	elseif i >relaxIter
		localArea = localArea - Dt * localAreaTimeStep;
		localArea = map.GD3.smoothDiffusionFFT(localArea, Dt, 10.0);
		localArea = localArea * InitialArea / map.surfaceIntegral(localArea); 
		protein = protein - Dt * proteinTimeStep;
		maskDepleted = protein < 0;
		protein(maskDepleted) = protein(maskDepleted) * 0.5;
		protein = map.GD3.smoothDiffusionFFT(protein, Dt, 10.0);
		protein = map.WENORK3Extend(protein, 100);
		%protein(protein<0) = 0;
		protein = protein * InitialArea * proteinCoverage / map.surfaceIntegral(protein); 
		%protein = protein - map.surfaceIntegral(protein)/InitialArea + proteinCoverage;
	end
	map.F = map.F - Dt * levelSetTimeStep;
	map.setDistance;

	localArea = map.WENORK3Extend(localArea, 100);
	protein = map.WENORK3Extend(protein, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bookkeeping
	ene_b = map.surfaceIntegral(0.5 * Kappa .* (MeanCurvature - SC).^2);
	ene_c = map.surfaceIntegral(0.5*Mu.*( (1./localArea - 1).^2 + (localArea - 1).^2));
	ene = ene_b + ene_c;

	array_t = [array_t time];
	array_b = [array_b; ene_b];
	array_c = [array_c; ene_c];

	fprintf('iter:%d, ene_b:%4.5f, ene_c:%4.5f, ar:%+4.5f, vol:%+4.5f, rd: %4.5f\n',i, ene_b, ene_c, DiffArea, DiffVolume, ReducedVolume)

	if i>1 && mod(i,1)==0
		clf(FIG)

		subplot(2,4,[7 8])
		area(array_t, [array_b array_c]);
		titleStr = [ sprintf('%5d: %.3e, ene_b:%.3f, ene_c:%.3f',i,time,ene_b,ene_c) ];
		title(titleStr)

		subplot(2,4,4)
		xslice = ceil(map.GD3.ncols / 2);
		F = reshape(map.F(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Y = reshape(map.GD3.Y(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Z = reshape(map.GD3.Z(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		contour(Y,Z,F,[0,0],'blue','LineWidth',3)
		axis equal
		titleStr = [sprintf(' shift:(%ld,%ld,%ld)', ...
				sign(x_shift),sign(y_shift),sign(z_shift))];
		title(titleStr)

		ax1 = subplot(2,4,[1 2 5 6]);
		titleStr = [ sprintf(' pC:%.3f, RTS:%.3f ', proteinCoverage,RelativeTimeScale) ];
		%map.plotField(0,localArea,0.0)
		map.plotField(0,protein,0.0)
		%map.plotField(0,MeanCurvature-SC,0.0)
		map.GD3.DrawBox
		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)
		%zoom(2.0)
		%set(ax1,'xlim',[-0.5,0.5],'ylim',[-0.5,0.5],'zlim',[-0.5,0.5])
		set(ax1,'xlim',[-0.80,0.80],'ylim',[-0.80,0.80],'zlim',[-0.80,0.80])

		%subplot(2,2,[1,3])
		ax3 = subplot(2,4,3);
		zoom reset
		titleStr = [ sprintf(' rd:%.3f, mu:%.3f ', ReducedVolume,Mu) ];
		map.plotField(0,localArea,0.0)
		%map.plotField(0,protein,0.0)
		%map.plotField(0,MeanCurvature-SC,0.0)
		map.GD3.DrawBox
		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)

		set(ax3,'xlim',[-0.80,0.80],'ylim',[-0.80,0.80],'zlim',[-0.80,0.80])
		%zoom(ax3,2.0)
		%zoom(ax1,2.0)
		%linkaxes([ax1,ax3],'xyz');
		%linkprop([ax1,ax3],'XLim','YLim','ZLim');

		drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if Archived 
			FIG.InvertHardcopy = 'off'; % preseve background color
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	end

	if mod(i,5)==0
		x_shift = 0; y_shift = 0; z_shift = 0;
		map.F = circshift(map.F, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		map.setDistance
		%if mod(i,20)==0
			map.F = map.WENO5RK3Reinitialization(map.F,200);
		%end

		localArea = circshift(localArea, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		localArea = map.WENORK3Extend(localArea,100);

		protein = circshift(protein, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		protein = map.WENORK3Extend(protein,100);
	end
end

%simu.simulationEnd
%SD.NE.processImage(30,'pear_protein_pinch')














