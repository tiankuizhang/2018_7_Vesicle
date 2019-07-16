% MeshFree Approximation Methods with Matlab Program30.2.
% PointCloud3D_PUCS

addpath(fullfile('PointCloud','kdtree','kdtree','lib'));
addpath(fullfile('WOBJ'))
addpath(fullfile('PointCloud'))

wf = @(e,r) r.^4.*(5*spones(r)-4*r);

rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
ep = 1;

npu = 32;
neval = 64;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xshift = (max(xe(:)) + min(xe(:))) / 2;
yshift = (max(ye(:)) + min(ye(:))) / 2;
zshift = (max(ze(:)) + min(ze(:))) / 2;
x = gpuArray(xe) - xshift;
y = gpuArray(ye) - yshift;
z = gpuArray(ze) - zshift;
F = gpuArray(reshape(Pf,neval,neval,neval));
Grid = SD.GD3(x,y,z);
map = SD.SDF3(Grid,x,y,z,F);

map.A = sqrt(x.^2 + y.^2 + (z-0.05).^2) - 0.03;

Ds = gather(map.GD3.Ds);
Dt = Ds.^2 * 1
for i = 1:100
	map.GPUsetCalculusToolBox
	levelSetTimeStep = - map.FRegularization(true);
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, 1.0);
	map.F = map.F - Dt * levelSetTimeStep;
%	map.F = map.GD3.smoothDiffusionFFT(map.F, Dt, .01);
end
for i = 1:200
	map.GPUsetCalculusToolBox
	levelSetTimeStep = - map.FRegularization(true);
	levelSetTimeStep = map.GD3.smoothFFT(levelSetTimeStep, Dt, 1.0);
	map.F = map.F - Dt * levelSetTimeStep;
	map.F = map.GD3.smoothDiffusionFFT(map.F, Dt, .01);
end

%map.F = map.ENORK2Reinitialization(map.F, 100);
%map.F = map.WENORK3Reinitialization(map.F, 100);
map.setDistance
map.F = map.WENO5RK3Reinitialization(map.F, 100);
%map.plotSurface(0,1,'r','none'); view(0,90)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map.setDistance
map.GPUsetCalculusToolBox
map.GPUAsetCalculusToolBox

FIG = figure('Name','Bunny','Position',[10 10 1600 800]);
iso = linspace(-0.1, 0.1, 50);

ax1 = subplot(1,2,1);
map.plotIsoField(iso, map.A, true)
axis equal vis3d
set(gca,'Color','k')
view(0,90)
axisLim = [-0.10 0.10];
set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
title('before')
colorbar off


map.A = map.WENO5RK3Extend(map.A, 100);
%map.A = map.WENORK3Extend(map.A, 100);
%map.A = map.ENORK2ClosetPointSurfaceRedistance(map.A,100,50);
%map.A = map.WENO5RK3ClosetPointSurfaceRedistance(map.A,50,100);
map.A = map.WENORK3ClosetPointSurfaceRedistance(map.A,20,100);
ax2 = subplot(1,2,2);
map.plotIsoField(iso, map.A, true)
axis equal vis3d
set(gca,'Color','k')
view(0,90)
axisLim = [-0.10 0.10];
set(gca,'xlim',axisLim,'ylim',axisLim,'zlim',axisLim)
title('after')
colorbar off

set(ax1, 'Position', [-0.2 0 1 1])
set(ax2, 'Position', [0.2 0 1 1])

































