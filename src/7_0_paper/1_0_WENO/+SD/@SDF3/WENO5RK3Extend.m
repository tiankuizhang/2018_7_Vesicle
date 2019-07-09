% extend a scalar field off the interface % boundary node values are fixed with interpolation
% upwind scheme used to define normal and therefore extend velocity
% ENO scheme used to calculate forward and backward derivative
% new extension scheme with subcell boundary fix
% first calculate the intersection of surface with grid lines
% then use ENO quadratic interpolation to calculate values at those points
% ENO derivatives near boundary will modified accordingly

function NewC = WENO5RK3Extend(obj, C, iteration)

	NewC = C;

	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	[fx, fy, fz] = feval(obj.WENO5RK3_upwind_normal, ...
			fx, fy, fz, obj.F, obj.xpr, obj.xpl, obj.ypf, obj.ypb, obj.zpu, obj.zpd, ...
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

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.WENO5RK3_boundary_interpolant, ...
			cpr,cpl,cpf,cpb,cpu,cpd,obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	deltat = 0.3 * obj.min_dist;

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.WENO5RK3_extend_step,step,deltat,NewC,vx,vy,vz,...
				obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp1 = NewC - step;

		step = feval(obj.WENO5RK3_extend_step,step,deltat,Ctmp1,vx,vy,vz,...
				obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;

		step = feval(obj.WENO5RK3_extend_step,step,deltat,Ctmp0_5,vx,vy,vz,...
				obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp1_5 = Ctmp0_5 - step;

		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
	end
	
end
