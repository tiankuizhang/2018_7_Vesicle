% calculate where the zero level set surface cut grid lines, i.e.
% set values for xpr etc
function setDistance(obj)

	obj.xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	obj.xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	obj.ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	obj.ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	obj.zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	obj.zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

	[obj.xpr, obj.xpl, obj.ypf, obj.ypb, obj.zpu, obj.zpd] = feval(obj.sixth_boundary_correction, ...
			obj.xpr, obj.xpl, obj.ypf, obj.ypb, obj.zpu, obj.zpd, obj.F, obj.GD3.NumElt, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

	minx = min(obj.xpr,obj.xpl);
	miny = min(obj.ypf,obj.ypb);
	minz = min(obj.zpu,obj.zpd);

	obj.min_dist = min(minx, min(miny,minz));

	xrange = obj.GD3.X(obj.xpr < obj.GD3.Dx | obj.xpl < obj.GD3.Dx);
	yrange = obj.GD3.Y(obj.ypf < obj.GD3.Dy | obj.ypb < obj.GD3.Dy);
	zrange = obj.GD3.Z(obj.zpu < obj.GD3.Dz | obj.zpd < obj.GD3.Dz);

	xmin = min(xrange); xmax = max(xrange);
	ymin = min(yrange); ymax = max(yrange);
	zmin = min(zrange); zmax = max(zrange);

	obj.Box = [xmin,xmax,ymin,ymax,zmin,zmax];

end
