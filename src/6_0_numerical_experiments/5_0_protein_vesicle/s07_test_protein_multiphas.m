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





















