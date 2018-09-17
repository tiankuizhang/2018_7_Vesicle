% test the extend scheme 

load('Extend.mat')

grid = SD.GD3(x,y,z);
map = SD.SDF3(grid,x,y,z,F);

map.GPUsetCalculusToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
clf

subplot(1,2,1)
%map.plot
map.plotField(0,map.MeanCurvature)

subplot(1,2,2)
tmp = map.WENORK3Extend(map.MeanCurvature,100);
map.plotField(0,tmp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Newton's method gives out of range value
% ind: 10076882 (66,65,62)
% ind: 1089599 (63,65,67)

% big difference between cubic and ENO
% ind:1007551 (63,64,62)

%s = map.GD3.Dx

%v0 = F(66,63,62);
%v1 = F(66,64,62);
%v2 = F(66,65,62);
%v3 = F(66,66,62);

%v0 = F(63,64,67);
%v1 = F(63,65,67);
%v2 = F(63,66,67);
%v3 = F(63,67,67);

%ind1 = 10076882;
%	
%	obj = map;
%	Fgpu = obj.F;
%
%	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
%	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
%	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
%	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
%	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
%	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
%
%	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.cubic_boundary_correction, ...
%			xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, obj.GD3.NumElt, ...
%		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
%			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
%
%	xpr1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
%	xpl1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
%	ypf1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
%	ypb1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
%	zpu1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
%	zpd1 = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
%	[xpr1, xpl1, ypf1, ypb1, zpu1, zpd1] = feval(obj.ENORK2_boundary_correction, ...
%			xpr1, xpl1, ypf1, ypb1, zpu1, zpd1, Fgpu, obj.GD3.NumElt, ...
%		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
%			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
%
%
%
%	diff = xpr - xpr1;
%	maxd = max(abs(diff(:)));
%	[ind1,ind2,ind3] = ind2sub(obj.GD3.Size, find(abs(diff)==maxd));
%	xpr(ind1,ind2,ind3)
%	xpr1(ind1,ind2,ind3)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
