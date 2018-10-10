% MeshFree Approximation Methods with Matlab Program30.1.
% PointCloud2D
% Scripts that fits a curve to 2D point cloud

rbf = @(e,r) exp(-(e*r).^2); ep = 3.5;
N = 81;
neval = 40;
t = 2*pi*haltonseq(N,1); dsites = [cos(t) sin(t)];
x = (2+sin(t)).*cos(t); y = (2+cos(t)).*sin(t);
nx = (2+cos(t)).*cos(t) - sin(t).^2;
ny = (2+sin(t)).*sin(t)-cos(t).^2;
dsites = [x y]; normals = [nx ny];

bmin = min(dsites,[],1); bmax = max(dsites,[],1);
bdim = max(bmax-bmin);

delta = bdim/100;

dsites(  N+1:2*N,:) = dsites(1:N,:) + delta*normals;
dsites(2*N+1:3*N,:) = dsites(1:N,:) - delta*normals;

rhs = [zeros(N,1); ones(N,1); -ones(N,1)];

ctrs = dsites;

bmin = min(dsites,[],1); bmax = max(dsites,[],1);

xgrid = linspace(bmin(1),bmax(1),neval);
ygrid = linspace(bmin(2),bmax(2),neval);
[xe,ye] = meshgrid(xgrid,ygrid);
epoints = [xe(:) ye(:)];
DM_eval = DistanceMatrix(epoints,ctrs);
EM = rbf(ep,DM_eval);
DM_data = DistanceMatrix(dsites,ctrs);
IM = rbf(ep,DM_data);
Pf = EM * (IM\rhs);

figure; hold on; view([-30,30])
plot3(dsites(:,1),dsites(:,2),rhs,'r.','markersize',20);
mesh(xe,ye,reshape(Pf,neval,neval));
axis tight; hold off

figure;hold on
plot(dsites(1:N,1),dsites(1:N,2),'bo');
contour(xe,ye,reshape(Pf,neval,neval),[0 0],'r');
hold off
