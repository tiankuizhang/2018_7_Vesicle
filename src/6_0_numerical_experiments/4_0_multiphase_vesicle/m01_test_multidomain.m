% test scheme to initialize a multidomain multiphase vesicle

radius = 0.98; ra = 1.5; xmax = radius*ra; xmin = -xmax;
GridSize = [64,64,64];

rd = 0.98;
alpha = pi/10;
domain = [...
			0,		pi/2,		alpha;...
			0,		-pi/2,		alpha;...
			0,		0,			alpha;...
			pi/2,	0,		alpha;...
			pi,		0,			alpha;...
			-pi/2,	0,			alpha;...
			];

[x,y,z,F,A,volume] = SD.Shape.MultiDomainSphere([xmin,xmax],GridSize,radius,rd,domain);
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);
map.A = A;


map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F,100);
map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);

map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox

subplot(1,2,1)
map.plotField(0,map.AHeaviside,0.01)
colorbar off
axis vis3d equal
map.GD3.DrawBox
set(gca,'Color','k')

subplot(1,2,2)
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
