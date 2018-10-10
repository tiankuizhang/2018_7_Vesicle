% MeshFree Approximation Methods with Matlab Program30.2.
% PointCloud3D_PUCS

wf = @(e,r) r.^4.*(5*spones(r)-4*r);

rbf = @(e,r) r.^4.*(5*spones(r)-4*r);
ep = 1;

npu = 8;

neval = 25;

% load position and normal data from ply file
Data3D_Bunny3 = fullfile('..','Objects','ply','stanford_bunny','bunny','reconstruction',...
		'bun_zipper_res4_n.ply');
ptCloud = pcread(Data3D_Bunny3);
dsites = ptCloud.Location;
normals = ptCloud.Normal;
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

bmin = min(dsites,[],1); bmax = max(dsites,[],1);
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
clellradius = 1/wep;

DM_eval = DistanceMatrixCSRBF(epoints,cellctrs,wep);





