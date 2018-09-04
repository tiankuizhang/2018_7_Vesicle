% test ENO-RK2 reinitialization scheme

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
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2)-0.6;

F = fun(x, y, z);

map = SD.SDF3(grid, x, y, z, F);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

obj = map;

Fgpu = obj.F;

%AF = obj.ENORK2Extend(obj.GD3.Z,100);
AF = obj.GD3.Z-0.3;
OldAF = AF;

xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

% AF instead of Fgpu should be used!!!!
[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
		xpr, xpl, ypf, ypb, zpu, zpd, AF, obj.GD3.NumElt, ...
	    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

% calculate the normal to the surface upwindly
fx = zeros(obj.GD3.Size, 'gpuArray');
fy = zeros(obj.GD3.Size, 'gpuArray');
fz = zeros(obj.GD3.Size, 'gpuArray');

[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
		fx, fy, fz, Fgpu, xpr, xpl, ypf, ypb, zpu, zpd, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

fgradient = sqrt(fx.^2+fy.^2+fz.^2);
flat = fgradient < 1e-6;
fgradient(flat) = 1e-6;

nx = fx ./ fgradient;
ny = fy ./ fgradient;
nz = fz ./ fgradient;

% calculate the sign of the original auxiliary level set function
Sign = zeros(obj.GD3.Size, 'gpuArray');
Sign(AF>0) = 1.;
Sign(AF<0) = -1.;

%
minx = min(xpr,xpl);
miny = min(ypf,ypb);
minz = min(zpu,zpd);
deltat = 0.3 * min(minx, min(miny,minz));

step = zeros(obj.GD3.Size, 'gpuArray');

tic
for i=1:100
	step = feval(obj.ENORK2_surface_redistance_step,step,AF,Sign,deltat,nx,ny,nz,...
		   		xpr,xpl,ypf,ypb,zpu,zpd,...	
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	AFtmp = AF - step;
	step = feval(obj.ENORK2_surface_redistance_step,step,AFtmp,Sign,deltat,nx,ny,nz,...
		   		xpr,xpl,ypf,ypb,zpu,zpd,...	
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	AF = (AF + AFtmp - step) / 2;

end
toc

AF = obj.ENORK2Extend(AF,100);

%	step = feval(obj.ENORK2_surface_redistance_step,step,AF,Sign,deltat,nx,ny,nz,...
%		   		xpr,xpl,ypf,ypb,zpu,zpd,...	
%		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
%				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure
%
%subplot(1,2,1)
%obj.plotIsoField([0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7],OldAF)
%
%subplot(1,2,2)
%obj.plotIsoField([0.1,0.2,0.3,0.4,0.5,0.6,0.7,-0.1,-0.2,-0.3,-0.4,-0.5,-0.6,-0.7],AF)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experiment with plotting intersection of two implicit surfaces

iso = 0.3;
[faces,verts,colors] = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,0,obj.GD3.Z);
mask = colors>iso;

outcount = sum(mask(faces),2);
cross = (outcount == 2) | (outcount == 1);
cross_tris = faces(cross,:);

%obj.plotSurface(0,1,'Green',1)
%ct = patch('Vertices',verts,'Faces',cross_tris,'EdgeColor',[1 .5 0],'FaceColor',[.5 1 .5]);

out_vert = mask(cross_tris);
flip = sum(out_vert,2) == 1;
out_vert(flip,:) = 1-out_vert(flip,:);

ntri = size(out_vert,1);
overt = zeros(ntri,3);

for i=1:ntri
	v1i = find(~out_vert(i,:));
	v2i = 1 + mod(v1i,3);
	v3i = 1 + mod(v1i+1,3);
	overt(i,:) = cross_tris(i,[v1i v2i v3i]);
end

u = (iso - colors(overt(:,1))) ./ (colors(overt(:,2)) - colors(overt(:,1)));
v = (iso - colors(overt(:,1))) ./ (colors(overt(:,3)) - colors(overt(:,1)));

uverts = repmat((1-u),[1 3]).*verts(overt(:,1),:) + repmat(u,[1 3]).*verts(overt(:,2),:);
vverts = repmat((1-v),[1 3]).*verts(overt(:,1),:) + repmat(v,[1 3]).*verts(overt(:,3),:);

x = nan(3,ntri);
x(1,:) = uverts(:,1)';
x(2,:) = vverts(:,1)';
y = nan(3,ntri);
y(1,:) = uverts(:,2)';
y(2,:) = vverts(:,2)';
z = nan(3,ntri);
z(1,:) = uverts(:,3)';
z(2,:) = vverts(:,3)';

obj.plotSurface(0,1,'Green',1)
h = line(x(:),y(:),z(:),'Color','red','LineWidth',3);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





























