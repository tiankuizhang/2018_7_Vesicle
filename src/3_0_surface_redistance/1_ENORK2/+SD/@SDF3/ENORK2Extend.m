% extend a scalar field off the interface
% boundary node values are fixed with interpolation
% upwind scheme used to define normal and therefore extend velocity
% ENO scheme used to calculate forward and backward derivative

function NewC = ENORK2Extend(obj, C, iteration)

	NewC = C;

	% calculate distance to the neighborhood node without crossing the interface
	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, obj.F, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
			fx, fy, fz, obj.F, xpr, xpl, ypf, ypb, zpu, zpd, ...
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
	
	% boundary nodes
	Boundary = ( obj.F .* obj.F(obj.GD3.oxo) < 0 ) | ...
			   ( obj.F .* obj.F(obj.GD3.oXo) < 0 ) | ...
			   ( obj.F .* obj.F(obj.GD3.yoo) < 0 ) | ...
			   ( obj.F .* obj.F(obj.GD3.Yoo) < 0 ) | ...
			   ( obj.F .* obj.F(obj.GD3.ooz) < 0 ) | ...
			   ( obj.F .* obj.F(obj.GD3.ooZ) < 0 ) ; 
	
	% closest points on boundary from boundary nodes 
	x_shift = obj.GD3.X(Boundary) - nx(Boundary) .* obj.F(Boundary);
	y_shift = obj.GD3.Y(Boundary) - ny(Boundary) .* obj.F(Boundary);
	z_shift = obj.GD3.Z(Boundary) - nz(Boundary) .* obj.F(Boundary);
	
	NewC(Boundary) = interp3(obj.GD3.X, obj.GD3.Y, obj.GD3.Z, C, ...
			x_shift, y_shift, z_shift, 'linear');
	
	step = zeros(obj.GD3.Size, 'gpuArray');
	deltat = 0.3 * min(obj.GD3.Dx, min(obj.GD3.Dy,obj.GD3.Dz));
	
	for i=1:iteration
		step = feval(obj.ENORK2_extend_step, step, NewC, Boundary, vx, vy, vz, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp = NewC - deltat * step;
		step = feval(obj.ENORK2_extend_step, step, Ctmp, Boundary, vx, vy, vz, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		NewC = (NewC + Ctmp - deltat * step) / 2;
	end
	
end
