% reinitialization with WENO and RK3 scheme for iteration times
% see 2008_Chene Second-Order Accurate Computation of Curvatures in a Level Set Framework
function NewF = WENORK3Reinitialization(obj,F,iteration)
	
	Fgpu = F;
	
	mask = F<0;
	deltat = zeros(obj.GD3.Size, 'gpuArray');

	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.cubic_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.WENORK3_reinitialization_step, step, Fgpu, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ftmp = Fgpu - step;
		step = feval(obj.WENORK3_reinitialization_step, step, Ftmp, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Fgpu = (Fgpu + Ftmp - step) / 2;
	end

	%obj.F = Fgpu;
	NewF = Fgpu;

end
