% test new scheme to account for protein dependent properties for mutiphase vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters
iteration = 5000; 
relaxIterC = 5000;
relaxIterP = 5000;
GridSize = [64,64,64]; 
Kappa0 = 1.0; Kappa1 = 0.0; % bending modulus for Ld phase
Kappa0Lo = 1.0; Kappa1Lo = 0.0; % bending modulus for Lo phase
KappaL = 70; % isotropic line tension
KappaG = 0; % difference in Gaussian bending rigidity: Ld - Lo
C0 = 0; C1 = -1; proteinCoverage = 1.0;
Mu = 1000; % incompressibility of vesicle
MinimumTimeStep = 0.0;
RelativeTimeScale = 1; % relative drag coefficient for protein motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha1 = pi/9;
beta1 = -1;
alpha2 = pi/10;
alpha3 = pi/8;
beta3 = pi/4;
domain1 = [...
			0,		pi/2,		alpha1,beta1;...
			0,		-pi/2,		alpha1,beta1;...
			0,		0,			alpha1,beta1;...
			pi/2,	0,			alpha1,beta1;...
			pi,		0,			alpha1,beta1;...
			-pi/2,	0,			alpha1,beta1;...
			];
domain2 = [...
			pi/4,pi/4,	alpha2,-1;	3*pi/4,pi/4,	alpha2,-1; ...
			-pi/4,pi/4,	alpha2,-1;	-3*pi/4,pi/4,	alpha2,-1;...
			pi/4,-pi/4,	alpha2,-1;	3*pi/4,-pi/4,	alpha2,-1; ...
			-pi/4,-pi/4,alpha2,-1;	-3*pi/4,-pi/4,	alpha2,-1;...
			];
domain3 = [...
			0,		pi/2,		alpha3,beta3;...
			0,		-pi/2,		alpha3,beta3;...
			];
domain4 = [...
			pi/4,pi/4,	alpha2,alpha2;	3*pi/4,pi/4,	alpha2,alpha2; ...
			-pi/4,pi/4,	alpha2,alpha2;	-3*pi/4,pi/4,	alpha2,alpha2;...
			pi/4,-pi/4,	alpha2,alpha2;	3*pi/4,-pi/4,	alpha2,alpha2; ...
			-pi/4,-pi/4,alpha2,alpha2;	-3*pi/4,-pi/4,	alpha2,alpha2;...
			];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius = 0.98; ra =2.0; xmax = radius*ra; xmin = -xmax; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multidomain pinching without collision
rd = 0.80;
%domain = [0,pi/2,0.53,pi/4];
domain = domain4;
Pressure = - 200; ConsereVol = false;
%ConsereVol = true;
CFLNumber = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
[x,y,z,F,A,volume] = SD.Shape.MultiDomainSphere2([xmin,xmax],GridSize,radius,rd,domain);
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);
map.A = A;

map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox

InitialArea = 4*pi*radius^2*1.0;
InitialVolume = (4*pi/3)*radius^3;
expectedVolume = volume;
AreaNegative = map.AcalArea;
AreaPositive = InitialArea - AreaNegative;

fprintf('initial area: %4.5f, expected volume: %4.5f\n', InitialArea, expectedVolume)
% name and size of figure
FIG = figure('Name','MultiPhase Vesicle','Position',[10 10 1600 800])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamics
time = 0;
array_t = [];	% time array
array_ld = [];	% bending energy for ld phase
array_lo = [];	% bending energy for lo phase
array_c = [];	% area incompressibility
array_el = [];	% line energy
array_eg = [];	% Guassian bending energy

localArea = ones(map.GD3.Size,'gpuArray');
protein = zeros(map.GD3.Size,'gpuArray');
for i = 1:iteration
	map.GPUsetCalculusToolBox
	map.GPUAsetCalculusToolBox

	z_shift = - (map.Box(5) + map.Box(6));
	y_shift = - (map.Box(3) + map.Box(4));
	x_shift = - (map.Box(1) + map.Box(2));

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - expectedVolume) / expectedVolume;
	ReducedVolume = (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);
	CurrentNegativeArea = map.AcalArea;
	CurrentPositiveArea = CurrentArea - CurrentNegativeArea;
	DiffPhaseArea = 100 * (CurrentNegativeArea - AreaNegative) / AreaNegative; 

	KappaLd = Kappa0 + Kappa1 * protein;
	KappaLo = Kappa0Lo + Kappa1Lo * protein;
	Kappa = map.BivalueField(KappaLd, KappaLo); % bending modulus field
	SC = C0 + C1 .* protein; % spontaneous curvature field

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% now calculate normal and tangential forces
	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	GaussianCurvature = map.WENORK3Extend(map.GaussianCurvature,100);
	localTension = 0.5*Mu .* ... % due to incompressibility
		( (1.0./localArea - 1).*(3.0./localArea - 1) + (1 - localArea.^2) );
	proTension = - 0.5 * Kappa1 * protein .* (MeanCurvature - SC).^2 ...
					+ C1 * protein .* Kappa .* (MeanCurvature - SC);

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

	% line forces in N direction
	[~,MeanCurvatureDn] = map.AsurfaceDerivative(MeanCurvature);
	GeodesicTorsion = map.WENORK3Extend(map.GeodesicTorsion,100);
	[GeodesicTorsionDt,~] = map.AsurfaceDerivative(GeodesicTorsion);

	LineSpeedN =  KappaL * map.NormalCurvature ...
				- KappaG * GeodesicTorsionDt ...
				+ (KappaLd - KappaLo) .* MeanCurvatureDn; ...
	LineSpeedN = map.ENORK2Extend(LineSpeedN,100);

	% line forces in n direction
	GeodesicCurvature = map.GeodesicCurvature;
	LineSpeedn =  KappaL * GeodesicCurvature ...
				- KappaG * GaussianCurvature ...
				- 0.5 * (KappaLd - KappaLo) .* MeanCurvature.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% (minus) total force in N direction
	NormalSpeedBend = NormalBendSpeed ...
		- LineSpeedN .* map.ADiracDelta .* map.AGradMag; 

	% determine the appropriate time step
	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	%Dt = max(MinimumTimeStep,Dt);
	time = time + Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% solve for tension and pressure to constrain total area and volume
	c11 = InitialArea; 
	c12 = map.AsurfaceIntegral(MeanCurvature); c21 = c12;
	c13 = map.surfaceIntegral(MeanCurvature.*map.AHeaviside); c31 = c12 + c13;
	c23 = - map.calLength; 
	c22 = map.AsurfaceIntegral(MeanCurvature.^2) - c23;
	c33 = map.surfaceIntegral(MeanCurvature.^2 .* map.AHeaviside) ;
	c32 = c33 + c22 + c23;

	tmp = map.GD3.LimitField(expectedVolume - CurrentVolume, 0.001*InitialVolume*CFLNumber);
	s1 = tmp / Dt + map.surfaceIntegral(NormalSpeedBend);

	sLine = map.LineIntegral(LineSpeedn);
	tmp = map.GD3.LimitField(AreaNegative-CurrentNegativeArea,0.01*InitialArea*CFLNumber);
	s2 = - tmp / Dt + ...
		map.AsurfaceIntegral(MeanCurvature.*NormalSpeedBend) + sLine;

	tmp = map.GD3.LimitField(InitialArea-CurrentArea,0.01*InitialArea*CFLNumber);
	s3 = - tmp / Dt + map.surfaceIntegral(MeanCurvature.*NormalSpeedBend);


	if ConsereVol 
		TP = [c11,c12,c13;c21,c22,c23;c31,c32,c33] \ [s1;s2;s3];
		Pressure = TP(1);
		TensionNegative = TP(2);
		TensionPositive = TP(3) + TP(2);
		Tension = map.BivalueField(TensionNegative, TensionPositive);
	else
		TP = [c22,c23;c32,c33] \ [s2-Pressure*c21; s3-Pressure*c31];
		TensionNegative = TP(1);
		TensionPositive = TP(2) + TP(1);
		Tension = map.BivalueField(TensionNegative, TensionPositive);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the complete normal speed 
	normalSpeed = Tension .* MeanCurvature + Pressure - NormalSpeedBend;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for the level set function
	maxKappa = max(abs(Kappa(mask)));
	levelSetTimeStep = normalSpeed.*map.FGradMag - maxKappa * map.FRegularization(true);
	levelSetTimeStep = normalSpeed.*map.FGradMag;
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, maxKappa);
	levelSetTimeStep = map.ENORK2Extend(levelSetTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for the auxilary level set function
	AnormalSpeed = LineSpeedn + TensionPositive - TensionNegative;
	[rgNormal, rgSD] = map.ARegularization(false);
	DN = 1.0; Dn = 1.0;
	AnormalSpeed = AnormalSpeed - DN * rgNormal - Dn * rgSD;
	AnormalSpeed = map.GD3.smoothFFT(AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% divergence of the flow field
	normalSpeed = map.GD3.smoothFFT(normalSpeed, Dt, maxKappa);
	Divergence =  - MeanCurvature .* normalSpeed ...
				+ map.GD3.GradDotGrad(0.5*MeanCurvature.^2, Kappa) ...
				+ 0.5*MeanCurvature.^2 .* map.GD3.Laplacian(Kappa) ...
				- map.GD3.GradDotGrad(MeanCurvature, Kappa.*SC) ...
				- MeanCurvature .* map.GD3.Laplacian(Kappa.*SC) ...
				+ map.GD3.Laplacian(0.5*Kappa.*SC.^2 + localTension + proTension);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (minus) time rate of change for localArea, i.e., sum of numerical fluxs
	laFlux = localArea .* Divergence + map.GD3.WENODotGrad(tvx,tvy,tvz,localArea);
	maxTv = max(abs(tv(mask)));
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
	pmaxTv = max(abs(pTv(mask)));
	proteinTimeStep = map.GD3.smoothDiffusionFFT(proFlux, Dt, pmaxTv);
	proteinTimeStep = map.WENORK3Extend(proteinTimeStep, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time step the system
	if i == relaxIterC
		localArea = ones(map.GD3.Size,'gpuArray');
	elseif i >relaxIterC
		localArea = localArea - Dt * localAreaTimeStep;
		localArea = map.GD3.smoothDiffusionFFT(localArea, Dt, 10.0);
		localArea = localArea * InitialArea / map.surfaceIntegral(localArea); 
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time step the system
	if i == relaxIterP
		expectedVolume = InitialVolume * 0.90;
		protein = proteinCoverage * ones(map.GD3.Size,'gpuArray');
		MinimumTimeStep = 1 * 1e-6;
	elseif i >relaxIterP
		protein = protein - Dt * proteinTimeStep;
		maskDepleted = protein < 0;
		protein(maskDepleted) = protein(maskDepleted) * 0.1;
		protein = map.GD3.smoothDiffusionFFT(protein, Dt, 10.0);
		protein = map.WENORK3Extend(protein, 100);
		%protein(protein<0) = 0;
		protein = protein * InitialArea * proteinCoverage / map.surfaceIntegral(protein); 
		%protein = protein - map.surfaceIntegral(protein)/InitialArea + proteinCoverage;
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	map.F = map.F - Dt * levelSetTimeStep;
	map.setDistance
	map.A = map.A - Dt * AnormalSpeed;

	localArea = map.WENORK3Extend(localArea, 100);
	protein = map.WENORK3Extend(protein, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bookkeeping
	ene_ld = map.AsurfaceIntegral(0.5 * Kappa .* (MeanCurvature - SC).^2);
	ene_lo = map.surfaceIntegral(0.5*Kappa.*(MeanCurvature-SC).^2.*map.AHeaviside);
	ene_c = map.surfaceIntegral(0.5*Mu.*( (1./localArea - 1).^2 + (localArea - 1).^2));
	ene_l = - KappaL * c23;
	ene_g = KappaG * map.LineIntegral(map.GeodesicCurvature);

	ene = ene_ld + ene_lo + ene_c + ene_l + ene_g;

	array_t = [array_t time];
	array_ld = [array_ld; ene_ld];
	array_lo = [array_lo; ene_lo];
	array_c = [array_c; ene_c];
	array_el = [array_el; ene_l];
	array_eg = [array_eg; ene_g];

	fprintf('iter: %5d, ene: %4.5f, ar: %+4.5f, vol: %+4.5f, rd: %4.5f, pe: %+4.5f\n', ...
			i, ene, DiffArea, DiffVolume, ReducedVolume, DiffPhaseArea)

	if i == 1
		ene0 = gather(ene);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if i>1 && mod(i,5)==0
		clf(FIG)

		subplot(2,4,8)
		area(array_t, [array_ld array_lo array_c array_el array_eg array_c]);
		ylim([0 ene0*1.2])
		titleStr = [ sprintf('%5d: %.3e, ene_b:%.3f, ene_c:%.3f', ...
				i,time,ene_ld+ene_lo,ene_c) ];
		title(titleStr)

		subplot(2,4,4)
		xslice = ceil(map.GD3.ncols / 2);
		Fslice = reshape(map.F(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Aslice = reshape(map.A(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Y = reshape(map.GD3.Y(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Z = reshape(map.GD3.Z(:,xslice,:), [map.GD3.mrows,map.GD3.lshts]);
		Fpositive = Fslice; Fpositive(Aslice<0) = nan;
		Fnegative = Fslice; Fnegative(Aslice>0) = nan;
		contour(Y,Z,Fpositive,[0,0],'blue','LineWidth',3)
		hold on
		contour(Y,Z,Fnegative,[0,0],'red','LineWidth',3)
		hold off
		axis equal
		titlestr = [sprintf('shift:(%1d,%1d,%1d)',sign(x_shift),sign(y_shift),sign(z_shift))];
		title(titlestr)

		ax3 = subplot(2,4,3);
		%ax1 = subplot(2,4,[1 2 5 6]);
		%ax7 = subplot(2,4,7);
		titleStr = [ sprintf('P pC:%.3f, RTS:%.3f ', proteinCoverage,RelativeTimeScale) ];
		%map.plotField(0,localArea,0.0)
		map.plotField(0,protein,0.0)
		colormap(gca,[parula]); caxis([0 2])
		%map.plotField(0,MeanCurvature-SC,0.0)
		map.GD3.DrawBox
		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)

		%ax1 = subplot(2,4,[1 2 5 6]);
		%ax3 = subplot(2,4,3);
		ax7 = subplot(2,4,7);
		zoom reset
		titleStr = [ sprintf('LA rd:%.3f, mu:%.3f ', ReducedVolume,Mu) ];
		map.plotField(0,localArea,0.0)
		colormap(gca,[parula]); caxis([0.9 1.1])
		map.GD3.DrawBox
		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)

		ax1 = subplot(2,4,[1 2 5 6]);
		%ax3 = subplot(2,4,3);
		%ax7 = subplot(2,4,7);
		zoom reset
		titleStr = [ sprintf(' rd:%.3f, mu:%.3f ', ReducedVolume,Mu) ];
		mask = map.GD3.X > map.GD3.Y; 
		Field = map.AHeaviside; 
		%Field(mask) = nan;
		map.plotField(0,Field,0.0); 
		alpha(0.6)
		colormap(gca,[1,0,0;0,0,1]); colorbar off; 
		map.GD3.DrawBox
		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
		axis vis3d equal
		set(gca,'Color','k')
		title(titleStr)
		set(ax1,'xlim',[-1.2,1.2],'ylim',[-1.2,1.2],'zlim',[-1.2,1.2])

		drawnow

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if mod(i,5)==0
		%x_shift = 0; y_shift = 0; z_shift = 0;
		map.F = circshift(map.F, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		map.setDistance
		%if mod(i,20)==0
			map.F = map.WENO5RK3Reinitialization(map.F,200);
		%end

		map.A = circshift(map.A, [sign(y_shift),sign(x_shift),sign(z_shift)]);
%		map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);

		localArea = circshift(localArea, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		localArea = map.WENORK3Extend(localArea,100);

		protein = circshift(protein, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		protein = map.WENORK3Extend(protein,100);
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

























