% test new scheme to account for protein dependent properties for mutiphase vesicle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulation parameters
iteration = 5000; relaxIter = 1000;
GridSize = [64,64,64]; 
Kappa0 = 1.0; Kappa1 = 0.0; % bending modulus for Ld phase
Kappa0Lo = 1.0; Kappa1Lo = 0.0; % bending modulus for Lo phase
KappaL = 70; % isotropic line tension
KappaG = 0; % difference in Gaussian bending rigidity: Ld - Lo
C0 = 0; C1 = -1; proteinCoverage = 1.0;
Mu = 1000; % incompressibility of vesicle
CFLNumber = 0.2;
MinimumTimeStep = 0.0;
RelativeTimeScale = 0.1; % relative drag coefficient for protein motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
radius = 0.98; ra =2.0; xmax = radius*ra; xmin = -xmax; rd = 0.87;
domain = [0,pi/2,0.53,-pi/4];
Pressure = - 200; ConsereVol = false;
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
	Kappa = map.BivalueField(KappaLd, KappaLo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	
end




















