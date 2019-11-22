%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%simu = SD.Simulation(mfilename, 'bunny');
%simu.simulationStart
%pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MeshFree Approximation Methods with Matlab Program30.2.
% PointCloud3D_PUCS

addpath(fullfile('PointCloud','kdtree','kdtree','lib'));
addpath(fullfile('WOBJ'))
addpath(fullfile('PointCloud'))

wf = @(e,r) r.^4.*(5*spones(r)-4*r);

rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
ep = 1;

npu = 32;
neval = 96;

% load position and normal data from ply file
Data3D_Bunny3 = fullfile('Objects','ply','stanford_bunny','bunny','reconstruction',...
		'bun_zipper_res3_n.ply');
ptCloud = pcread(Data3D_Bunny3);
dsites = double(ptCloud.Location);
normals = double(ptCloud.Normal);

%OBJ = read_wobj();
%dsites = double(OBJ.vertices);
%normals = double(OBJ.vertices_normal);
N = size(dsites,1);

bmin = min(dsites,[],1); bmax = max(dsites,[],1);
bdim = max(bmax-bmin);
wep = npu/bdim;

withnormals = find(normals(:,1)|normals(:,2)|normals(:,3));
addpoints = length(withnormals);

delta = bdim/100;

dsites(N+1          :N+addpoints,:)   = dsites(withnormals,:) + delta*normals(withnormals,:);
dsites(N+addpoints+1:N+2*addpoints,:) = dsites(withnormals,:) - delta*normals(withnormals,:);

rhs = [zeros(N,1); ones(addpoints,1); -ones(addpoints,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% center grid
bmin = min(dsites,[],1); bmax = max(dsites,[],1);
bshift = (bmax + bmin) / 2;
N = size(dsites,1);
dsites = dsites - repmat(bshift, [N, 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bmin = min(dsites,[],1); bmax = max(dsites,[],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add buffer points 
ratio = 1.3;
bmin = ratio * bmin; bmax = ratio * bmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ctrs = dsites;

xgrid = linspace(bmin(1),bmax(1),neval);
ygrid = linspace(bmin(2),bmax(2),neval);
zgrid = linspace(bmin(3),bmax(3),neval);
[xe,ye,ze] = meshgrid(xgrid,ygrid,zgrid);
epoints = [xe(:) ye(:) ze(:)];

puxgrid = linspace(bmin(1),bmax(1),npu);
puygrid = linspace(bmin(2),bmax(2),npu);
puzgrid = linspace(bmin(3),bmax(3),npu);
[xpu,ypu,zpu] = meshgrid(puxgrid,puygrid,puzgrid);
cellctrs = [xpu(:) ypu(:) zpu(:)];
cellradius = 1/wep;

DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);
SEM = wf(wep,DM_eval);
SEM = spdiags(1./(SEM*ones(npu^3,1)),0,neval^3,neval^3)*SEM;

[tmp,tmp,datatree] = kdtree(dsites, []);
[tmp,tmp,evaltree] = kdtree(epoints,[]);
Pf = zeros(neval^3,1);
for j=1:npu^3
	[pts,dist,idx] = kdrangequery(datatree,cellctrs(j,:),cellradius);
	if (length(idx) > 0)
		DM_data = DistanceMatrixCSRBF(dsites(idx,:),ctrs(idx,:),ep);
		IM = rbf(ep,DM_data);
		[epts,edist,edix] = kdrangequery(evaltree,cellctrs(j,:),cellradius);
		DM_eval = DistanceMatrixCSRBF(epoints(edix,:),ctrs(idx,:),ep);
		EM = rbf(ep,DM_eval);
		localfit = EM * (IM\rhs(idx));
		Pf(edix) = Pf(edix) + localfit.*SEM(edix,j);
	end
end

%figure; hold on
%%plot3(dsites(1:N,1),dsites(1:N,2),dsites(1:N,3),'bo');
%pfit = patch(isosurface(xe,ye,ze,reshape(Pf,neval,neval,neval),0));
%isonormals(xe,ye,ze,reshape(Pf,neval,neval,neval),pfit);
%set(pfit,'FaceLighting','gouraud','FaceColor','red','EdgeColor','none');
%light('Position',[0 0 1],'Style','infinite');
%daspect([1 1 1]); view([0 90]);
%axis([bmin(1) bmax(1) bmin(2) bmax(2) bmin(3) bmax(3)]);
%%axis off; 
%hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu = SD.Simulation(mfilename, 'bunny');
simu.simulationStart
pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map.F = map.ENORK2Reinitialization(map.F, 100);

xshift = (max(xe(:)) + min(xe(:))) / 2;
yshift = (max(ye(:)) + min(ye(:))) / 2;
zshift = (max(ze(:)) + min(ze(:))) / 2;
x = gpuArray(xe) - xshift;
y = gpuArray(ye) - yshift;
z = gpuArray(ze) - zshift;
F = gpuArray(reshape(Pf,neval,neval,neval));
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);


Ds = gather(map.GD3.Ds);
Dt = Ds.^2 * 1
for i = 1:100
	map.GPUsetCalculusToolBox
	levelSetTimeStep = - 0.5 * map.FRegularization(true);
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, 1.0);
	map.F = map.F - Dt * levelSetTimeStep;
	map.F = map.GD3.smoothDiffusionFFT(map.F, Dt, .015);
	%if mod(i,10)==0
	%if true
	%	map.F = map.ENORK2Reinitialization(map.F, 10);
	%end
end
for i = 1:200
%	map.GPUsetCalculusToolBox
%	levelSetTimeStep = - map.FRegularization(true);
%	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, 1.0);
%	map.F = map.F - Dt * levelSetTimeStep;
	map.F = map.GD3.smoothDiffusionFFT(map.F, Dt, .01);
end
for i = 1:0
	map.GPUsetCalculusToolBox
	levelSetTimeStep = - map.FRegularization(true);
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, 1.0);
	map.F = map.F - Dt * levelSetTimeStep;
	map.F = map.GD3.smoothDiffusionFFT(map.F, Dt, .01);
	%if mod(i,10)==0
	%if true
	%	map.F = map.ENORK2Reinitialization(map.F, 10);
	%end
end

map.F = map.ENORK2Reinitialization(map.F, 100);
map.F = map.WENORK3Reinitialization(map.F, 100);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F, 100);
%map.plotSurface(0,1,'r','none'); view(0,90)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map.setDistance
map.GPUsetCalculusToolBox

%map.A = sqrt(x.^2 + y.^2 + (z-0.05).^2) - 0.03;
map.A = exp(x.^2 + y.^2 + (z-0.05).^2) - exp(0.03^2);
map.GPUAsetCalculusToolBox

FIG = figure('Name','Bunny','Position',[10 10 1600 800]);
%FIG = figure('Name','Bunny','Position',[10 10 800 400]);
iso = linspace(-0.1, 0.1, 50);

%ax1 = subplot(1,2,1);
%map.plotIsoField(iso, map.A, true)
%axis equal vis3d
%set(gca,'Color','k')
%view(0,90)
%axisLim = [-0.10 0.10];
%set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
%title('before')
%colorbar off


map.A = map.WENO5RK3Extend(map.A, 100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map.A = map.WENO5RK3ClosetPointSurfaceRedistance(map.A,200,100);
	obj = map;
	A = map.A;
	iter1 = 200;
	iter2 = 100;
	%% calculated for the reinitialization scheme
	% calculate distance to the neighborhood node without crossing the interface
	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.sixth_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, A, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));
	
	mask = A < 0;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%% calculated for the extend scheme
	% calculate distance to the neighborhood node without crossing the interface
	fxpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	fxpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	fypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	fypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	fzpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	fzpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	[fxpr, fxpl, fypf, fypb, fzpu, fzpd] = feval(obj.sixth_boundary_correction, ...
			fxpr, fxpl, fypf, fypb, fzpu, fzpd, obj.F, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	minx = min(fxpr,fxpl);
	miny = min(fypf,fypb);
	minz = min(fzpu,fzpd);
	fdeltat = 0.3 * min(minx, min(miny,minz));

	%fdeltat = fdeltat .* obj.F ./ sqrt(obj.F.^2 + obj.GD3.Dx);
	
	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	[fx, fy, fz] = feval(obj.WENO5RK3_upwind_normal, ...
			fx, fy, fz, obj.F, fxpr, fxpl, fypf, fypb, fzpu, fzpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	fgradient = max(sqrt(fx.^2+fy.^2+fz.^2), 1e-14);
	
	nx = fx ./ fgradient;
	ny = fy ./ fgradient;
	nz = fz ./ fgradient;

%	% tmp code
%	r = sqrt(obj.GD3.X.^2 + obj.GD3.Y.^2 + obj.GD3.Z.^2);
%	nx = obj.GD3.X./r;
%	ny = obj.GD3.Y./r;
%	nz = obj.GD3.Z./r;
	
	Sign = sign(obj.F);
	vx = Sign .* nx;
	vy = Sign .* ny;
	vz = Sign .* nz;
	
	% interpolate values for crossing points
	cpr = zeros(obj.GD3.Size, 'gpuArray');
	cpl = zeros(obj.GD3.Size, 'gpuArray');
	cpf = zeros(obj.GD3.Size, 'gpuArray');
	cpb = zeros(obj.GD3.Size, 'gpuArray');
	cpu = zeros(obj.GD3.Size, 'gpuArray');
	cpd = zeros(obj.GD3.Size, 'gpuArray');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	step = zeros(obj.GD3.Size, 'gpuArray');
	NewA = A;
	for i=1:iter1
		NewA = ENORK2Extend(obj,NewA,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

		step = feval(obj.WENO5RK3_reinitialization_step, step, NewA, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
		Atmp1 = NewA - step;
		Atmp1 = ENORK2Extend(obj,Atmp1,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);
	
		step = feval(obj.WENO5RK3_reinitialization_step, step, Atmp1, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
		Atmp2 = Atmp1 - step;

		Atmp0_5 = 3./4. * NewA + 1./4. * Atmp2;
		Atmp0_5 = ENORK2Extend(obj,Atmp0_5,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

		step = feval(obj.WENO5RK3_reinitialization_step, step, Atmp0_5, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
		Atmp1_5 = Atmp0_5 - step;

		NewA = 1./3. * NewA + 2./3. * Atmp1_5;

		%if mod(i,5) == 0
		if true
			map.A = NewA;
			clf
			map.plotIsoField(iso, map.A, true)
			axis equal vis3d
			set(gca,'Color','k')
			view(0,90)
			axisLim = [-0.10 0.10];
			set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
			title([ sprintf('iter: %3d', i) ])
			colorbar off

			set(gca, 'Position', [0 0 1 1])
			drawnow
			FIG.InvertHardcopy = 'off'; % preserve background color
			saveas(FIG, fullfile(simu.JPG, [sprintf('%05d',i),'isosurface','.jpg']))
			
		end

	end
	NewA = ENORK2Extend(obj,NewA,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

	map.A = NewA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simu.simulationEnd
simu.processImage(20)

%ax2 = subplot(1,2,2);
%map.plotIsoField(iso, map.A, true)
%axis equal vis3d
%set(gca,'Color','k')
%view(0,90)
%axisLim = [-0.10 0.10];
%set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
%title('after')
%colorbar off
%
%set(ax1, 'Position', [-0.2 0 1 1])
%set(ax2, 'Position', [0.2 0 1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NewC = ENORK2Extend(obj, C, iteration, xpr,xpl,ypf,ypb,zpu,zpd,vx,vy,vz,deltat,...
					cpr,cpl,cpf,cpb,cpu,cpd)

	NewC = C;

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.WENO5RK3_boundary_interpolant, ...
			cpr,cpl,cpf,cpb,cpu,cpd,xpr,xpl,ypf,ypb,zpu,zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
	
	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.WENO5RK3_extend_step,step,deltat,NewC,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
		Ctmp1 = NewC - step;

		step = feval(obj.WENO5RK3_extend_step,step,deltat,Ctmp1,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);

		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;

		step = feval(obj.WENO5RK3_extend_step,step,deltat,Ctmp0_5,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);

		Ctmp1_5 = Ctmp0_5 - step;

		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
	end
	
end
































