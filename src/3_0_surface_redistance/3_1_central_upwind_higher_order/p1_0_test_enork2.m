% create a 3D grid
xv = linspace(-1,1,64);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

% create a SDF3 instance
Radius = 0.6;
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-Radius;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);
map.GPUsetCalculusToolBox;

A = z - 0.3;
Extend = map.ENORK2Extend(A,100);
Re = map.ENORK2Extend( map.ENORK2Reinitialization(Extend,100), 100);
Sur = map.ENORK2Extend( map.ENORK2CentralUpwindSurfaceRedistance(Extend,100), 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the surface redistancing equation with the closet point method

obj = map;
C = Extend;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculated for the reinitialization scheme
% calculate distance to the neighborhood node without crossing the interface
xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
		xpr, xpl, ypf, ypb, zpu, zpd, Extend, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

minx = min(xpr,xpl);
miny = min(ypf,ypb);
minz = min(zpu,zpd);
deltat = 0.3 * min(minx, min(miny,minz));

mask = Extend < 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculated for the extend scheme
% calculate distance to the neighborhood node without crossing the interface
fxpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
fxpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
fypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
fypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
fzpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
fzpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

[fxpr, fxpl, fypf, fypb, fzpu, fzpd] = feval(obj.ENORK2_boundary_correction, ...
		fxpr, fxpl, fypf, fypb, fzpu, fzpd, obj.F, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

minx = min(fxpr,fxpl);
miny = min(fypf,fypb);
minz = min(fzpu,fzpd);

fdeltat = 0.3 * min(minx, min(miny,minz));

% calculate the extend velocity upwindly
fx = zeros(obj.GD3.Size, 'gpuArray');
fy = zeros(obj.GD3.Size, 'gpuArray');
fz = zeros(obj.GD3.Size, 'gpuArray');

[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
		fx, fy, fz, obj.F, fxpr, fxpl, fypf, fypb, fzpu, fzpd, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

fgradient = sqrt(fx.^2+fy.^2+fz.^2);
flat = fgradient < 1e-6;
fgradient(flat) = 1e-6;

nx = fx ./ fgradient;
ny = fy ./ fgradient;
nz = fz ./ fgradient;

vx = sign(obj.F) .* nx;
vy = sign(obj.F) .* ny;
vz = sign(obj.F) .* nz;

% interpolate values for crossing points
cpr = zeros(obj.GD3.Size, 'gpuArray');
cpl = zeros(obj.GD3.Size, 'gpuArray');
cpf = zeros(obj.GD3.Size, 'gpuArray');
cpb = zeros(obj.GD3.Size, 'gpuArray');
cpu = zeros(obj.GD3.Size, 'gpuArray');
cpd = zeros(obj.GD3.Size, 'gpuArray');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
step = zeros(obj.GD3.Size, 'gpuArray');
for i=1:100
	step = feval(obj.ENORK2_reinitiliaztion_step, step, C, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp = C - step;
	Ctmp = ENORK2Extend(obj,Ctmp,100,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,...
			cpr,cpl,cpf,cpb,cpu,cpd);

	step = feval(obj.ENORK2_reinitiliaztion_step, step, Ctmp, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp2 = Ctmp - step;
	Ctmp2 = ENORK2Extend(obj,Ctmp2,100,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,...
			cpr,cpl,cpf,cpb,cpu,cpd);

	C = (C + Ctmp2) / 2;

end
toc
Re = obj.ENORK2Extend(C, 100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = Extend;
tic
step = zeros(obj.GD3.Size, 'gpuArray');
for i=1:10
	step = feval(obj.ENORK2_reinitiliaztion_step, step, C, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp = C - step;
	Ctmp = ENORK2Extend(obj,Ctmp,5,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,...
			cpr,cpl,cpf,cpb,cpu,cpd);

	step = feval(obj.ENORK2_reinitiliaztion_step, step, Ctmp, mask, deltat, ...
			xpr, xpl, ypf, ypb, zpu, zpd, ...
	    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	Ctmp2 = Ctmp - step;
	Ctmp2 = ENORK2Extend(obj,Ctmp2,5,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,...
			 cpr,cpl,cpf,cpb,cpu,cpd);

	C = (C + Ctmp2) / 2;

end
toc
Sur = obj.ENORK2Extend(C,100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%map.A = Extend;
%map.AsetCalculusToolBox
%
%figure
%
%subplot(2,2,1)
%map.plotField(0,map.GeodesicCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,2)
%map.plotField(0,map.NormalCurvature)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,3)
%map.plotField(0,map.GeodesicTorsion)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);
%
%subplot(2,2,4)
%map.plotField(0,map.BPerpendicular)
%map.plotIsoField([-0.2,-0.1,0.1,0.2],map.A,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare effects of extendsion,reinitialization and surface redistance scheme on
% calculation of curvature

figure

iso = -0.5:0.1:0.5;

map.A = Re;
map.AsetCalculusToolBox

subplot(2,2,1)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);
subplot(2,2,2)
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);

map.A = Sur;
map.AsetCalculusToolBox
subplot(2,2,3)
map.plotField(0,map.ADiracDelta)
map.plotIsoField(iso,map.A,false);
subplot(2,2,4)
map.plotField(0,map.GeodesicCurvature)
map.plotIsoField(iso,map.A,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function NewC = ENORK2Extend(obj, C, iteration, xpr,xpl,ypf,ypb,zpu,zpd,vx,vy,vz,deltat,...
					cpr,cpl,cpf,cpb,cpu,cpd)

	NewC = C;

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.ENORK2_boundary_interpolate, ...
			cpr,cpl,cpf,cpb,cpu,cpd,xpr,xpl,ypf,ypb,zpu,zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.ENORK2_extend_step,step,deltat,NewC,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp = NewC - step;
		step = feval(obj.ENORK2_extend_step,step,deltat,Ctmp,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		NewC = (NewC + Ctmp - step) / 2;
	end
	
end








