% test new approach to regularize A field

iter = 80;

radius = 0.98; ra =2.0; xmax = radius*ra; xmin = -xmax; GridSize = [64,64,64];
rd = 0.87;
alpha1 = pi/9;
beta1 = -1;
alpha2 = pi/10;
alpha3 = pi/6;
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

%domain = [domain1;domain2];
%domain = domain2;
%domain = domain3;
%domain = [0,pi/2,pi/5,-pi/4];
domain = [0,pi/2,0.53,-pi/4]; % rd = 0.87

%domain = domain1;
%domain = [...
%			0,		pi/2,		alpha1,beta1;...
%			0,		-pi/2,		alpha1,beta1;...
%			0,		0,			alpha1*1.1,beta1;...
%			pi/2,	0,			alpha1,beta1;...
%			pi,		0,			alpha1,beta1;...
%			-pi/2,	0,			alpha1,beta1;...
%			];

Pressure = - 200; ConsereVol = false;
%ConsereVol = true;

[x,y,z,F,A,volume] = SD.Shape.MultiDomainSphere2([xmin,xmax],GridSize,radius,rd,domain);
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);
map.A = A;

% bending mudulus
KappaB = 1.0; % bending rigidity for Ld phase
KappaBLo = 1 * KappaB; % bending rigidity for Lo phase
KappaL = 70; % isotropic line tension
KappaG = 0; % difference in Gaussian bending rigidity: Ld - Lo


map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,100,50);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox
InitialArea = 4*pi*radius^2*1.0;
%InitialArea = 4*pi*radius^2;
InitialVolume = (4*pi/3)*radius^3;
%expectedVolume = InitialVolume;
expectedVolume = volume;
AreaNegative = map.AcalArea;
AreaPositive = InitialArea - AreaNegative;

fprintf('initial area: %4.5f, expected volume: %4.5f\n', InitialArea, expectedVolume)
% name and size of figure
FIG = figure('Name','MultiPhase Vesicle','Position',[10 10 1600 800])


% dynamics
time = 0;
array_t = [];
array_ld = [];
array_lo = [];
array_el = [];
array_eg = [];

CFLNumber = 1;
localArea = ones(map.GD3.Size,'gpuArray');
filterWidth = gather(map.GD3.Ds)*5.0;
for i = 1:iter
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% surface forces
	KappaBField = map.BivalueField(KappaB, KappaBLo); 
	MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);
	GaussianCurvature = map.WENORK3Extend(map.GaussianCurvature,100);
	%MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(map.MeanCurvature); 
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(MeanCurvature); 
	NormalBendSpeed = KappaBField .* ( MeanCurvatureSurfaceLaplacian + ...
	 0.5 * MeanCurvature.^3 - 2 * MeanCurvature .* GaussianCurvature );
	NormalBendSpeed = map.ENORK2Extend(NormalBendSpeed, 100);
	%mask = abs(map.F)<2*map.GD3.Ds;
	%MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	% line forces in N direction
	[~,MeanCurvatureDn] = map.AsurfaceDerivative(MeanCurvature);
	GeodesicTorsion = map.WENORK3Extend(map.GeodesicTorsion,100);
	[GeodesicTorsionDt,~] = map.AsurfaceDerivative(GeodesicTorsion);

	LineSpeedN =  KappaL * map.NormalCurvature ...
				- KappaG * GeodesicTorsionDt ...
				+ (KappaB - KappaBLo) * MeanCurvatureDn; ...
	LineSpeedN = map.ENORK2Extend(LineSpeedN,100);

	% line forces in n direction
	%GeodesicCurvature = map.AENORK2Extend(map.GeodesicCurvature,50,100,50);
	%GeodesicCurvature = map.WENORK3Extend(map.GeodesicCurvature,100);
	GeodesicCurvature = map.GeodesicCurvature;
	LineSpeedn =  KappaL * GeodesicCurvature ...
				- KappaG * GaussianCurvature ...
				- 0.5 * (KappaB - KappaBLo) * MeanCurvature.^2;
	%LineSpeedn = map.AENORK2Extend(LineSpeedn, 50, 100, 50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% (minus) total force in N direction
	NormalSpeedBend = NormalBendSpeed ...
		- LineSpeedN .* map.ADiracDelta .* map.AGradMag; 

	% determine the appropriate time step
	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% calculate local Lagrange multiplier
	normalSpeedBendSmoothed = map.GD3.smoothFFT(NormalSpeedBend, Dt, ...
			0.5*max(KappaB,KappaBLo) );
	localArea = map.WENORK3Extend(localArea,100);
	totalArea = map.surfaceIntegral(localArea);
	localArea = localArea * CurrentArea / totalArea; % rescale to conserve total area
	Alpha = map.surfaceIntegral(MeanCurvature.^2) / InitialArea;
	%rhs = ((localArea - 1)./localArea)/Dt ;
	%rhs = ((localArea - 1)./localArea)/Dt - MeanCurvature.*NormalBendSpeed;
	%rhs = ((localArea - 1)./localArea)/Dt - MeanCurvature.*NormalSpeedBend;
	%rhs = ((localArea - 1)./localArea)/Dt - MeanCurvature.*normalSpeedBendSmoothed;
	% if 0.01 is changed to 0.05,01 cylindrical symmetry will be broken
	%rhs = map.GD3.LimitField( ((localArea - 1)./localArea)/Dt, 0.01/Dt) ...
	%	- MeanCurvature.*NormalBendSpeed;
	%rhs = map.GD3.LimitField( ((localArea - 1)./localArea)/Dt, 0.01/Dt) ...
	%	- MeanCurvature.*NormalSpeedBend;
	%rhs = map.GD3.LimitField( ((localArea - 1)./localArea)/Dt, 0.1/Dt) ...
	%	- MeanCurvature.*normalSpeedBendSmoothed;
	%rhs = - MeanCurvature.*NormalBendSpeed;
	%[localTension,residual] = localLagrangeMultiplier(map,rhs,Dt,Alpha);
	%localTension = map.WENORK3Extend(localTension,100);
	%residual = map.WENORK3Extend(residual,100);
	localTension = zeros(map.GD3.Size,'gpuArray');

	NormalSpeedBend = NormalSpeedBend - localTension .* MeanCurvature;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% solve for tension and pressure to constrain total area and volume
	c11 = InitialArea; 
	%c11 = AreaNegative; 
	c12 = map.AsurfaceIntegral(MeanCurvature); c21 = c12;
	c13 = map.surfaceIntegral(MeanCurvature.*map.AHeaviside); c31 = c12 + c13;
	c23 = - map.calLength; 
	c22 = map.AsurfaceIntegral(MeanCurvature.^2) - c23;
	c33 = map.surfaceIntegral(MeanCurvature.^2 .* map.AHeaviside) ;
	c32 = c33 + c22 + c23;

	tmp = map.GD3.LimitField(expectedVolume - CurrentVolume, 0.001*InitialVolume);
	s1 = tmp / Dt + map.surfaceIntegral(NormalSpeedBend);

	sLine = map.LineIntegral(LineSpeedn);
	%tmp = map.GD3.LimitField(AreaNegative-CurrentNegativeArea,0.0001*CurrentNegativeArea);
	tmp = map.GD3.LimitField(AreaNegative-CurrentNegativeArea,0.01*InitialArea);
	s2 = - tmp / Dt + ...
		map.AsurfaceIntegral(MeanCurvature.*NormalSpeedBend) + sLine;

	tmp = map.GD3.LimitField(InitialArea-CurrentArea,0.01*InitialArea);
	s3 = - tmp / Dt + map.surfaceIntegral(MeanCurvature.*NormalSpeedBend);

%	TP = [c11,c12,c13;c21,c22,c23;c31,c32,c33] \ [s1;s2;s3];
%	PressureNegative = TP(1);
%	TensionNegative = TP(2);
%	TensionPositive = TP(3) + TP(2);
%	Tension = map.BivalueField(TensionNegative, TensionPositive);
%	Pressure = map.BivalueField(PressureNegative, 0);

%	if DiffVolume < 1 
%		ConsereVol = true;
%		KappaL = 70;
%	end
	
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
%
	if i == 2
		[c11,c12,c13;c21,c22,c23;c31,c32,c33] 
		[s1;s2;s3]
		keyboard
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	% now calculate normal Speed for both level sets
	normalSpeed = Tension .* MeanCurvature - NormalSpeedBend + Pressure;
	% regularize both speed fields
	normalSpeedSmoothed = map.GD3.smoothFFT( normalSpeed.*map.FGradMag, Dt, ...
			0.5*max(KappaB,KappaBLo) );
	% it is import to use KappaBLo instead of KappaB 
	%normalSpeedSmoothed = smoothGMRES(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);

	AnormalSpeed = LineSpeedn + TensionPositive - TensionNegative;
	%AnormalSpeed = map.GD3.smoothFFT(AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
	%AnormalSpeed = smoothDiffusionFFT(map, AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
	%AnormalSpeed = smoothGMRES(map, AnormalSpeed.*map.AGradMag, Dt, 0.5);
	AnormalSpeed = map.GD3.smoothFFT(AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
	%AnormalSpeed = map.AENORK2Extend(AnormalSpeed,50,100,50); % reduce asymmteric error

	[rgNormal, rgSD] = map.ARegularization(false);
	DN = 1.0;
	Dn = 1.0;

	AnormalSpeed = AnormalSpeed - DN * rgNormal - Dn * rgSD;
	AnormalSpeed = map.GD3.smoothFFT(AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
	%AnormalSpeed = AnormalSpeed - map.GD3.smoothFFT(DN * rgNormal + Dn * rgSD, ...
	%		Dt, max(DN, Dn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% update c field  given normalSpeed and AnormalSpeed
	%Divergence = residual - (localArea - 1) ./ (localArea * Dt) ...
	%	- Tension .* MeanCurvature.^2 - Pressure .* MeanCurvature;
	% change due to divergence of velocity field
	%normalSpeed0 = (Tension + localTension) .* MeanCurvature + Pressure ...
	%	- NormalBendSpeed; % ignore line forcw while compute divergence
	%Divergence = - MeanCurvature .* normalSpeed0 + map.GD3.Laplacian(localTension);
	Divergence = - MeanCurvature .* normalSpeedSmoothed + map.GD3.Laplacian(localTension);
	%localArea = ones(map.GD3.Size,'gpuArray');
	%Divergence = residual - Tension .* MeanCurvature.^2 - Pressure .* MeanCurvature;
	%Divergence = - Tension .* MeanCurvature.^2 - Pressure .* MeanCurvature;
	%localArea = localArea - Dt * Divergence.*localArea;
	%localArea = map.WENORK3Extend(localArea,100);
	% advection velocity and advection term
	[vx,vy,vz] = map.GD3.Gradient(localTension);
	[cx,cy,cz] = map.GD3.Gradient(localArea);
	Advection = map.GD3.DotProduct(vx,vy,vz,cx,cy,cz);
	%mask = abs(map.A)<4*map.GD3.Ds;
	%Advection = zeros(map.GD3.Size, 'gpuArray');
	%Advection = feval(map.advection_step,Advection,vx,vy,vz,localArea,...
	%		map.GD3.mrows,map.GD3.ncols,map.GD3.lshts,...
	%		map.GD3.Dx,map.GD3.Dy,map.GD3.Dz);
	%Advection(mask) = 0;
	DcDt = map.WENORK3Extend(Divergence.*localArea+Advection,100);
	%DcDt = smoothDiffusionFFT(map,DcDt,Dt,0.5);
	%DcDt = map.ENORK2Extend(DcDt, 100);
	localArea = localArea - Dt * DcDt;
	%localArea = imgaussfilt3(localArea, filterWidth); % fileter to remove high order error
	%localArea = map.WENORK3Extend(localArea,100);
	%localArea = smoothDiffusionFFT(map,localArea,Dt,10);
	localArea = map.GD3.smoothDiffusionFFT(localArea,Dt,1);
	%localArea = map.WENORK3Extend(localArea,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	

	% timestep level set function
	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance
	map.A = map.A - Dt * AnormalSpeed;
	%map.A = map.WENORK3Extend(map.A, 50);
	%map.A = map.WENORK3Reinitialization(map.A, 100);
	%map.A = map.WENORK3Extend(map.A, 50);
	localArea = map.WENORK3Extend(localArea,100);

	% do some bookkeeping
	ene_ld = 0.5 * KappaB * (c22 + c23);
	ene_lo = 0.5 * KappaBLo * c33;
	ene_l = - KappaL * c23;
	ene_g = KappaG * map.LineIntegral(map.GeodesicCurvature);
	ene = ene_ld + ene_lo + ene_l + ene_g;

	array_t = [array_t time];
	array_ld = [array_ld; ene_ld];
	array_lo = [array_lo; ene_lo];
	array_el = [array_el; ene_l];
	array_eg = [array_eg; ene_g];

	fprintf('iter: %5d, ene: %4.5f, ar: %+4.5f, vol: %+4.5f, rd: %4.5f, pe: %+4.5f\n', ...
			i, ene, DiffArea, DiffVolume, ReducedVolume, DiffPhaseArea)
	

	%if mod(i,3)==0 || i==2
	if  i>1
		timeStr = [sprintf('%04d: %0.5e, %0.5f', i,time,ene)];

		clf(FIG)

		subplot(2,2,[1,3])
		titlestr = [ sprintf(' rd:%.3f, kd:%.1f,ko:%.1f,kl:%.1f,kg:%.1f ', ...
				ReducedVolume,KappaB,KappaBLo,KappaL,KappaG) ];
		title(titlestr)
		Field = map.AHeaviside; %Field(map.GD3.X>2*map.GD3.Ds)=nan;
		map.plotField(0,Field,0.0);
		colormap(gca,[1,0,0;0,0,1]);
		%map.plotField(0,map.AHeaviside,0.0);colormap(gca,[1,0,0;0,0,1]);
		colorbar off; 
		%map.plotSurface(0,1,'Green','black');textZ = gather(map.GD3.zmin);
		%map.plotField(0,normalSpeedSmoothed,0.5)
		%map.plotField(0,map.AHeaviside,0.0);colormap([1,0,0;0,0,1]);
		%colorbar off; camlight; lighting gouraud
		%map.plotField(0,map.ADiracDelta,0.01)
		%map.plotField(0,Tension,0.01)
		%map.plotField(0,localTension+Tension,0.01)
		%map.plotField(0,localTension,0.01)
		%map.plotField(0,residual,0.01)
		%map.plotField(0,Divergence,0.01)
		map.GD3.DrawBox

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		set(gca,'Color','k')

		zoom(1.0)

		subplot(2,2,4)
		area(array_t, [array_ld array_lo array_el array_eg])
		titleStr = [ sprintf('%5d: %.3e, energy:%.3f',i,time,ene) ];
		title(titleStr)
		
%		subplot(2,2,2)
%		map.plotField(0,localArea,0.01)
%		map.GD3.DrawBox
%		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
%		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
%		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])
%		axis vis3d equal
%		set(gca,'Color','k')
%		titlestr = [sprintf('shift:(%1d,%1d,%1d)',sign(x_shift),sign(y_shift),sign(z_shift))];
%		title(titlestr)

		subplot(2,2,2)
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
		%extraStr = [sprintf('\n %.3e , %.3e',volumeChangelt, volumeChangelt*Dt)];
		%title([titlestr,extraStr])
		title(titlestr)

		drawnow


		if false 
			FIG.InvertHardcopy = 'off'; % preseve background color
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,5)==0

		map.F = circshift(map.F, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		map.setDistance
		map.F = map.WENO5RK3Reinitialization(map.F,200);

		map.A = circshift(map.A, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		%map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
		%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,20,30);

		% we need to shift every field!!!!
		localArea = circshift(localArea, [sign(y_shift),sign(x_shift),sign(z_shift)]);
		localArea = map.WENORK3Extend(localArea,100);
	end

	if mod(i,10)==0
		%map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,10,0);
	end

end





