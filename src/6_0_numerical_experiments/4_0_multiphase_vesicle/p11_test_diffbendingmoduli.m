% simulate a two phase vesicle conserving volume and area of different phases
% incomprssibility is not imposed right now
% a new scheme to impose constraints used
% we shall add GaussianBending energy in this script
% we shall try to match our results with Baumgart's experiment
% start from a sphere, keep phase area and total area fixed
% relax volume
% add different bending moduli


% rd: reduced volume, ration of area for the lipid disordered phase
%rd = 0.76; AreaRatioLd = 0.56; ra = 2;
rd = 0.84; AreaRatioLd = 0.18; ra = 2;
%rd = 0.98; AreaRatioLd = 0.89; ra = 1.5;
%rd = 0.98; AreaRatioLd = 0.95;

GridSize = [64,64,64];
radius = 0.98;
xmax = radius * ra; xmin = - xmax;

% bending mudulus
KappaB = 1.0; % bending rigidity for Ld phase
KappaBLo = 5.0; % bending rigidity for Lo phase
KappaL = 70; % isotropic line tension
KappaG = 0; % difference in Gaussian bending rigidity


iter = 1500;

[x,y,z,F,A,volume] = SD.Shape.BiphaseSphere([xmin,xmax],GridSize,radius,rd,AreaRatioLd);
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);
map.A = A;

map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,100,50);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox
InitialArea = 4*pi*radius^2;
InitialVolume = (4*pi/3)*radius^3;
expectedVolume = volume;
AreaNegative = InitialArea * AreaRatioLd;
AreaPositive = InitialArea * (1 - AreaRatioLd);


fprintf('initial area: %4.5f, expected volume: %4.5f\n', InitialArea, expectedVolume)
% name and size of figure
FIG = figure('Name','MultiPhase Vesicle','Position',[10 10 1600 800])

iso = gather(-20*map.GD3.Ds:4*map.GD3.Ds:20*map.GD3.Ds);

% position to show iteration number and accumulated time
textX = gather(map.GD3.xmin);
%textY = gather( (map.GD3.ymax + map.GD3.ymin)/2 );
textY = gather( map.GD3.ymin );
textZ = gather(map.GD3.zmin);

CFLNumber = 1;
filterWidth = gather(map.GD3.Ds)*5.0;

% dynamics
time = 0;
array_t = [];
array_eb = [];
array_el = [];
array_eg = [];

z_cen = map.surfaceIntegral(map.GD3.Z);
z_shift = - floor(z_cen / map.GD3.Dz);

VolumeReached = false;
KappaBField = KappaB;

for i = 1:iter
	%KappaL = i;
	map.GPUsetCalculusToolBox
	map.GPUAsetCalculusToolBox

	CurrentArea = map.calArea;
	DiffArea = 100 * (CurrentArea - InitialArea)/InitialArea;
	CurrentVolume = map.calVolume;
	DiffVolume = 100 * (CurrentVolume - expectedVolume) / expectedVolume;
	ReducedVolume = (3*CurrentVolume/4/pi) * (4*pi/CurrentArea)^(3/2);
	CurrentNegativeArea = map.AcalArea;
	CurrentPositiveArea = CurrentArea - CurrentNegativeArea;
	DiffPhaseArea = 100 * (CurrentNegativeArea - AreaNegative) / AreaNegative; 

	if abs(DiffVolume) < 0.1
		VolumeReached = true;
	end
	if VolumeReached
		KappaBField = map.BivalueField(KappaB, KappaBLo); 
	end


	% extend mean curvature away from surface
	map.MeanCurvature = map.WENORK3Extend(map.MeanCurvature,100);

	% surface Laplacian of mean curvature and numerical Hamiltonian
	MeanCurvatureSurfaceLaplacian = map.GD3.Laplacian(map.MeanCurvature); 
	NormalSpeedBend = KappaBField .* (MeanCurvatureSurfaceLaplacian + ...
	 0.5 * map.MeanCurvature .* (map.MeanCurvature.^2 - 4 * map.GaussianCurvature) );
	NormalSpeedBend = map.ENORK2Extend(NormalSpeedBend, 100);

	%mask = abs(map.F)<2*map.GD3.Ds;
	%MaxSpeedBend = max(abs(NormalSpeedBend(mask)));

	NormalCurvature = map.ENORK2Extend(map.NormalCurvature, 100);
	LineNormalSpeed = KappaL * map.GeodesicCurvature - KappaG * map.GaussianCurvature;
	LineNormalSpeed = map.AENORK2Extend(LineNormalSpeed, 50, 100, 50);

	NormalSpeedBend = NormalSpeedBend ...
	 - KappaL * NormalCurvature .* map.ADiracDelta .* map.AGradMag ;

	mask = abs(map.F)<2*map.GD3.Ds;
	MaxSpeedBend = max(abs(NormalSpeedBend(mask)));
	%MaxSpeedBend = max( max(abs(NormalSpeedBend(mask))), 10*MaxSpeedBend );

	Dt = CFLNumber * map.GD3.Ds / MaxSpeedBend;
	time = time + Dt;

	% now solve for tension and pressure to enfore total area and volume
	c11 = InitialArea; 
	c12 = map.AsurfaceIntegral(map.MeanCurvature); c21 = c12;
	c13 = map.surfaceIntegral(map.MeanCurvature.*map.AHeaviside); c31 = c13;
	c23 = - map.calLength; c32 = c23;
	c22 = map.AsurfaceIntegral(map.MeanCurvature.^2) - c23;
	c33 = map.surfaceIntegral(map.MeanCurvature.^2 .* map.AHeaviside) - c23;

	tmp = expectedVolume - CurrentVolume;
	tmp = sign(tmp) * min(abs(tmp), 0.001*CurrentVolume);
	s1 = tmp / Dt + map.surfaceIntegral(NormalSpeedBend);
	%s1 = (expectedVolume - CurrentVolume) / Dt + ...
	%	map.surfaceIntegral(NormalSpeedBend);

	sLine = map.LineIntegral(LineNormalSpeed);
	s2 = - (AreaNegative - CurrentNegativeArea) / Dt + ...
		map.AsurfaceIntegral(map.MeanCurvature.*NormalSpeedBend) + sLine;
	s3 = - (AreaPositive - CurrentPositiveArea) / Dt + ...
		map.surfaceIntegral(map.MeanCurvature.*NormalSpeedBend.*map.AHeaviside)-sLine;...

	TP = [c11,c12,c13;c21,c22,c23;c31,c32,c33] \ [s1;s2;s3];
	Pressure = TP(1);
	TensionNegative = TP(2);
	TensionPositive = TP(3);
	Tension = map.BivalueField(TensionNegative, TensionPositive);

	% now calculate normal Speed
	normalSpeed = Tension .* map.MeanCurvature - NormalSpeedBend + Pressure;

	% time step level set function
	normalSpeedSmoothed = smoothFFT(map, normalSpeed.*map.FGradMag, Dt, 0.5*KappaBLo);
	% it is import to use KappaBLo instead of KappaB 
	%normalSpeedSmoothed = smoothGMRES(map, normalSpeed.*map.FGradMag, Dt, 0.5);
	normalSpeedSmoothed = map.ENORK2Extend(normalSpeedSmoothed, 100);

	%% time step the auxilary level set function
	AnormalSpeed = LineNormalSpeed + TensionPositive - TensionNegative;

	%AnormalSpeed = smoothFFT(map, AnormalSpeed.*map.AGradMag, Dt, 0.5);
	AnormalSpeed = smoothDiffusionFFT(map, AnormalSpeed.*map.AGradMag, Dt, 0.5*KappaL);
	%AnormalSpeed = smoothGMRES(map, AnormalSpeed.*map.AGradMag, Dt, 0.5);

	% timestep level set function
	map.F = map.F - Dt * normalSpeedSmoothed;
	map.setDistance
	map.A = map.A - Dt * AnormalSpeed;
	%map.A = map.WENORK3Extend(map.A, 50);
	%map.A = map.WENORK3Reinitialization(map.A, 100);
	map.A = map.WENORK3Extend(map.A, 50);

	ene_b = KappaB * (c22 + c33 + 2*c23);
	if VolumeReached
		ene_b = KappaB * (c22 + c23) + KappaBLo * ( c33 + c23);
	end
	
	ene_l = - KappaL * c23;
	ene_g = KappaG * map.LineIntegral(map.GeodesicCurvature);
	ene = ene_b + ene_l + ene_g;

	array_t = [array_t time];
	array_eb = [array_eb; ene_b];
	array_el = [array_el; ene_l];
	array_eg = [array_eg; ene_g];

	fprintf('iter: %5d, ene: %4.5f, ar: %+4.5f, vol: %+4.5f, rd: %4.5f, pe: %+4.5f\n', ...
			i, ene, DiffArea, DiffVolume, ReducedVolume, DiffPhaseArea)
	

	if mod(i,20)==0 || i==2
		timeStr = [sprintf('%04d: %0.5e, %0.5f', i,time,ene)];

		clf(FIG)

		subplot(2,2,[1,3])
		titlestr = [ sprintf( ' rd:%.2f,kl:%.2f,gr:(%d,%d,%d),\n zcen:%.2f, zshift: %d', rd,KappaL,GridSize(1),GridSize(2),GridSize(3),z_cen,z_shift ) ];
		title(titlestr)
		%map.plotSurface(0,1,'Green','black');textZ = gather(map.GD3.zmin);
		%map.plotField(0,normalSpeedSmoothed,0.5)
		%map.plotField(0,map.AHeaviside,0.01)
		map.plotField(0,Tension,0.01)
		map.GD3.DrawBox

		xticks([map.GD3.BOX(1),0,map.GD3.BOX(2)])
		yticks([map.GD3.BOX(3),0,map.GD3.BOX(4)])
		zticks([map.GD3.BOX(5),0,map.GD3.BOX(6)])

		axis vis3d equal
		th=text(textX, textY, textZ, timeStr, 'Color', 'y', 'FontSize', 14);
		set(th,'BackgroundColor', 'k', 'EdgeColor', 'w')

		zoom(1.0)

		subplot(2,2,4)
		area(array_t, [array_eb array_el array_eg])
		
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

		drawnow


		if false 
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
		end
	end

	if mod(i,10)==0

		z_cen = map.surfaceIntegral(map.GD3.Z);
		z_shift = - floor(z_cen / map.GD3.Dz);
	%	z_shift = sign(z_shift);

	%	map.F = circshift(map.F, [0,0,sign(z_shift)]);
	%	map.setDistance
		map.F = map.WENO5RK3Reinitialization(map.F,200);
		%map.GPUsetCalculusToolBox

	%	map.A = circshift(map.A, [0,0,sign(z_shift)]);
		map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
		%map.GPUAsetCalculusToolBox

		%map.F = map.ENORK2Reinitialization(map.F,100);
		%map.F = map.WENO5RK3Reinitialization(map.F,200);
		%map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,20,30);
		%map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
		
	end

end

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



