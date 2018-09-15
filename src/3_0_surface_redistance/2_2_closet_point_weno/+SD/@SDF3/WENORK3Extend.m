% extend a scalar field off the interface % boundary node values are fixed with interpolation
% upwind scheme used to define normal and therefore extend velocity
% ENO scheme used to calculate forward and backward derivative
% new extension scheme with subcell boundary fix
% first calculate the intersection of surface with grid lines
% then use ENO quadratic interpolation to calculate values at those points
% ENO derivatives near boundary will modified accordingly

function NewC = WENORK3Extend(obj, C, iteration)

	NewC = C;

	% calculate distance to the neighborhood node without crossing the interface
	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	%[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.cubic_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, obj.F, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	%[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
	[fx, fy, fz] = feval(obj.WENORK3_upwind_normal, ...
			fx, fy, fz, obj.F, xpr, xpl, ypf, ypb, zpu, zpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	fgradient = max(sqrt(fx.^2+fy.^2+fz.^2), 1e-6);
	
	nx = fx ./ fgradient;
	ny = fy ./ fgradient;
	nz = fz ./ fgradient;
	
	vx = sign(obj.F) .* nx;
	vy = sign(obj.F) .* ny;
	vz = sign(obj.F) .* nz;
	
	% interpolate values for crossing points
	cpr = zeros(obj.GD3.Size, 'gpuArray');
	cpl = zeros(obj.GD3.Size, 'gpuArray');
	cpf = zeros(obj.GD3.Size, 'gpuArray');
	cpb = zeros(obj.GD3.Size, 'gpuArray');
	cpu = zeros(obj.GD3.Size, 'gpuArray');
	cpd = zeros(obj.GD3.Size, 'gpuArray');

	%[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.ENORK2_boundary_interpolate, ...
	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.WENORK3_boundary_interpolate, ...
			cpr,cpl,cpf,cpb,cpu,cpd,xpr,xpl,ypf,ypb,zpu,zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		%step = feval(obj.WENORK3_extend_step,step,deltat,NewC,vx,vy,vz,...
		step = feval(obj.ENORK2_extend_step,step,deltat,NewC,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp1 = NewC - step;

		%step = feval(obj.WENORK3_extend_step,step,deltat,Ctmp1,vx,vy,vz,...
		step = feval(obj.ENORK2_extend_step,step,deltat,Ctmp1,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;

		%step = feval(obj.WENORK3_extend_step,step,deltat,Ctmp0_5,vx,vy,vz,...
		step = feval(obj.ENORK2_extend_step,step,deltat,Ctmp0_5,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp1_5 = Ctmp0_5 - step;

		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
	end
	
end
