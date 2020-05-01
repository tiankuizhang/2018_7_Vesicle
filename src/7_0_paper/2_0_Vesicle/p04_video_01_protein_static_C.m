%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('pear.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu = SD.Simulation(mfilename, 'protein_static_C');
simu.simulationStart
Archived = true;
pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters
iteration = 520;
GridSize = [64,64,64]; ReducedVolume0 = 0.65; VesicleTYPE = "o"; ratio = 0.2;
[x,y,z,~] = SD.Shape.Ellipsoid(GridSize,ReducedVolume0,VesicleTYPE,ratio);
Kappa0 = 1.0; Kappa1 = 0.0; C0 = -0.0; C1 = -1; proteinCoverage = 4.75;
%Kappa0 = 1.0; Kappa1 = 10.0; C0 = -0.0; C1 = 0; proteinCoverage = 1;
%Kappa0 = 1.0; Kappa1 = 10.0; C0 = -0.0; C1 = -4.5; proteinCoverage = 1;
Mu = 0; % incompressibility of vesicle
CFLNumber = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.GPUsetCalculusToolBox
InitialArea = map.calArea;
InitialVolume = map.calVolume;
InitialReducedVolume = (3*InitialVolume/4/pi) * (4*pi/InitialArea)^(3/2);
expectedVolume = InitialVolume; % no water crossing membrane

fprintf('initial area: %.3e, initial volume: %.3e, rd: %.3e\n', InitialArea, InitialVolume, InitialReducedVolume)
FIG = figure('Name','Vesicle','Position',[10 10 1600 800])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iteration
time = 0;
array_t = []; % time array
array_b = []; % bending energy
array_c = []; % area compressibility

localArea = ones(map.GD3.Size,'gpuArray');
protein = proteinCoverage * ones(map.GD3.Size,'gpuArray');
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
	% (minus) bending forces
	NormalBendSpeed = map.GD3.Laplacian( Kappa.*(MeanCurvature-SC) ) + ...
		+ Kappa .* (0.5*MeanCurvature.^3 - 2*MeanCurvature.*GaussianCurvature ...
					 + 2*SC.*GaussianCurvature - 0.5*SC.^2.*MeanCurvature) ...
		- MeanCurvature .* localTension ;
	NormalBendSpeed = map.ENORK2Extend(NormalBendSpeed,100);
	% tangential forces
	[tempx1,tempy1,tempz1] = map.GD3.Gradient(Kappa);
	[tempx2,tempy2,tempz2] = map.GD3.Gradient(Kappa.*SC);
	[tempx3,tempy3,tempz3] = map.GD3.Gradient(0.5*Kappa.*SC.^2 + localTension);
	tvx = 0.5 * MeanCurvature.^2 .* tempx1 - MeanCurvature .* tempx2 + tempx3;
	tvy = 0.5 * MeanCurvature.^2 .* tempy1 - MeanCurvature .* tempy2 + tempy3;
	tvz = 0.5 * MeanCurvature.^2 .* tempz1 - MeanCurvature .* tempz2 + tempz3;
	tv = sqrt(tvx.^2 + tvy.^2 + tvz.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select time step
	mask = abs(map.F) < 2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalBendSpeed(mask)));
	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contrain total area and volume
	c11 = map.surfaceIntegral(MeanCurvature.^2);
	c12 = map.surfaceIntegral(MeanCurvature); c21 = c12;
	c22 = CurrentArea;

	areaChangeRate = (InitialArea - CurrentArea) / Dt;
	volumeChangeRate = (InitialVolume - CurrentVolume) / Dt;
	s1 = - areaChangeRate + map.surfaceIntegral(NormalBendSpeed.*MeanCurvature);
	s2 = volumeChangeRate + map.surfaceIntegral(NormalBendSpeed);

	TP = [c11,c12;c21,c22] \ [s1;s2];
	Tension = TP(1); Pressure= TP(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the complete normal speed 
	normalSpeed = Tension .* MeanCurvature + Pressure - NormalBendSpeed;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for the level set function
	maxKappa = max(abs(Kappa(mask)));
	levelSetTimeStep = map.GD3.smoothFFT(normalSpeed.*map.FGradMag, Dt, maxKappa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divergence of the flow field
	%normalSpeedSmoothed = map.GD3.smoothFFT(normalSpeed, Dt, maxKappa);
	%Divergence =  - MeanCurvature .* normalSpeedSmoothed ...
	Divergence =  - MeanCurvature .* normalSpeed ...
				+ map.GD3.GradDotGrad(0.5*MeanCurvature.^2, Kappa) ...
				+ 0.5*MeanCurvature.^2 .* map.GD3.Laplacian(Kappa) ...
				- map.GD3.GradDotGrad(MeanCurvature, Kappa.*SC) ...
				- MeanCurvature .* map.GD3.Laplacian(Kappa.*SC) ...
				+ map.GD3.Laplacian(0.5*Kappa.*SC.^2 + localTension);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for localArea, i.e., sum of numerical fluxs
	laFlux = localArea .* Divergence + map.GD3.DotGrad(tvx,tvy,tvz,localArea);
	maxTv = max(abs(tv(mask)));
	localAreaTimeStep = map.GD3.smoothDiffusionFFT(laFlux, Dt, maxTv);
	localAreaTimeStep = map.WENORK3Extend(localAreaTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% protein tension
%	proTension = - 0.5 * Kappa1 * protein .* (MeanCurvature - SC).^2 ...
%					+ C1 * protein .* Kappa .* (MeanCurvature - SC);
%% divergence of protein motion field
%	proDivergence = Divergence ...
%		+ map.GD3.Laplacian(proTension) - MeanCurvature.^2 .* proTension;
%% additional tangential motion
%	[ptvx,ptvy,ptvz] = map.GD3.Gradient(proTension);
%% (minus) time rate of change for protein field
%	proFlux = protein .* proDivergence ...
%		+ map.GD3.DotGrad(tvx+ptvx, tvy+ptvy, tvz+ptvz, protein);
%	proteinTimeStep = map.GD3.smoothDiffusionFFT(proFlux, Dt, maxTv);
%	proteinTimeStep = map.WENORK3Extend(proteinTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% protein tension
	proTension = - 0.5 * Kappa1 * protein .* (MeanCurvature - SC).^2 ...
					+ C1 * protein .* Kappa .* (MeanCurvature - SC);
% divergence of protein motion field, only tangential effects are accounted for
	proDivergence = map.GD3.Laplacian(proTension);
% additional tangential motion
	[ptvx,ptvy,ptvz] = map.GD3.Gradient(proTension);
% (minus) time rate of change for protein field
	proFlux = protein .* proDivergence ...
		+ map.GD3.DotGrad(ptvx, ptvy, ptvz, protein);
	proteinTimeStep = map.GD3.smoothDiffusionFFT(proFlux, Dt, maxTv);
	proteinTimeStep = map.WENORK3Extend(proteinTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time step the system
	localArea = localArea - Dt * localAreaTimeStep;
	localArea = localArea * InitialArea / map.surfaceIntegral(localArea); 
	protein = protein - Dt * proteinTimeStep;
	protein = map.GD3.smoothDiffusionFFT(protein, Dt, 1.0);
	protein(protein<0) = 0;
	protein = protein * InitialArea * proteinCoverage / map.surfaceIntegral(protein); 
	%protein = protein - map.surfaceIntegral(protein)/InitialArea + proteinCoverage;
%	map.F = map.F - Dt * levelSetTimeStep;
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

	if i>1 & mod(i,5)==0
		clf(FIG)

		subplot(2,2,4)
		area(array_t, [array_b array_c]);
		titleStr = [ sprintf('%5d: %.3e, ene_b:%.3f, ene_c:%.3f',i,time,ene_b,ene_c) ];
		title(titleStr)

		%ax2 = subplot(2,2,2);
		ax2 = subplot(2,2,[1,3]);
		titleStr = [ sprintf(' rd:%.3f, mu:%.3f ', ...
			ReducedVolume,Mu) ];
		%map.plotField(0,localArea,0.0)
		map.plotField(0,protein,0.0)
		%map.plotField(0,protein,0.0)
		%map.plotField(0,MeanCurvature-SC,0.0)
		map.GD3.DrawBox
		caxis([0 17]) % for spontanous curvature
		%caxis([0 2.5]) % for bending moduli

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)
		set(ax2,'xlim',[-0.80,0.80],'ylim',[-0.80,0.80],'zlim',[-0.80,0.80])

		%ax13 = subplot(2,2,[1,3]);
		ax13 = subplot(2,2,2);
		titleStr = [ sprintf(' rd:%.3f, mu:%.3f ', ...
			ReducedVolume,Mu) ];
		%map.plotField(0,localArea,0.0)
		%map.plotField(0,protein,0.0)
		%map.plotField(0,protein,0.0)
		%map.plotField(0,MeanCurvature-SC,0.0)
		map.plotField(0,-MeanCurvature,0.0)
		map.GD3.DrawBox
		caxis([0 17])

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)
		set(ax13,'xlim',[-0.80,0.80],'ylim',[-0.80,0.80],'zlim',[-0.80,0.80])


		drawnow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		if Archived
			FIG.InvertHardcopy = 'off'; % preseve background color
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	pause(1)
	end

	if mod(i,5)==0
		map.F = circshift(map.F, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		map.setDistance
		%map.F = map.WENO5RK3Reinitialization(map.F,200);

		localArea = circshift(localArea, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		localArea = map.WENORK3Extend(localArea,100);

		protein = circshift(protein, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		protein = map.WENORK3Extend(protein,100);
	end


end

simu.simulationEnd
SD.NE.processImage(30,'protein_static_C')

