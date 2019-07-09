function NewC = ANENORK2Extend(obj, C, iter1, iter2, iter3)

	% extend C from surface F
	NewC = obj.ENORK2Extend(C, iter1);
	%NewC = C;

	% calculate distance to the neighborhood node without crossing the interface
	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, obj.A, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	% calculate the extend velocity upwindly
	Ax = zeros(obj.GD3.Size, 'gpuArray');
	Ay = zeros(obj.GD3.Size, 'gpuArray');
	Az = zeros(obj.GD3.Size, 'gpuArray');
	
	[Ax, Ay, Az] = feval(obj.ENORK2_upwind_normal, ...
			Ax, Ay, Az, obj.A, xpr, xpl, ypf, ypb, zpu, zpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	Agradient = sqrt(Ax.^2+Ay.^2+Az.^2);
	flat = Agradient < 1e-6;
	Agradient(flat) = 1e-6;
	
	nx = Ax ./ Agradient;
	ny = Ay ./ Agradient;
	nz = Az ./ Agradient;
	
	vx = sign(obj.A) .* nx;
	vy = sign(obj.A) .* ny;
	vz = sign(obj.A) .* nz;
	
	% interpolate values for crossing points
	cpr = zeros(obj.GD3.Size, 'gpuArray');
	cpl = zeros(obj.GD3.Size, 'gpuArray');
	cpf = zeros(obj.GD3.Size, 'gpuArray');
	cpb = zeros(obj.GD3.Size, 'gpuArray');
	cpu = zeros(obj.GD3.Size, 'gpuArray');
	cpd = zeros(obj.GD3.Size, 'gpuArray');

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.ENORK2_boundary_interpolate, ...
			cpr,cpl,cpf,cpb,cpu,cpd,xpr,xpl,ypf,ypb,zpu,zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
	
	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iter2
		step = feval(obj.ENORK2_extend_step,step,deltat,NewC,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp = NewC - step;
		step = feval(obj.ENORK2_extend_step,step,deltat,Ctmp,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		NewC = (NewC + Ctmp - step) / 2;
	end
	
	NewC = obj.ENORK2Extend(NewC, iter3);

end
