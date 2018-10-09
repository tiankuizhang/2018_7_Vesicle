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
