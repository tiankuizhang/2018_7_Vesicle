% redistancing A with the colset point methodology
% it seems that iter1 = 10, iter2 = 5 is good enough for our purposes
% A should be extended before being passed to this scheme
function NewA = WENORK3ClosetPointSurfaceRedistance(obj,A,iter1,iter2)

	%% calculated for the reinitialization scheme
	% calculate distance to the neighborhood node without crossing the interface
	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	%[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.cubic_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, A, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));
	
	mask = A < 0;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	%% calculated for the extend scheme
	% calculate distance to the neighborhood node without crossing the interface
	fxpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	fxpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	fypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	fypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	fzpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	fzpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	%[fxpr, fxpl, fypf, fypb, fzpu, fzpd] = feval(obj.ENORK2_boundary_correction, ...
	[fxpr, fxpl, fypf, fypb, fzpu, fzpd] = feval(obj.cubic_boundary_correction, ...
			fxpr, fxpl, fypf, fypb, fzpu, fzpd, obj.F, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	
	
	minx = min(fxpr,fxpl);
	miny = min(fypf,fypb);
	minz = min(fzpu,fzpd);
	fdeltat = 0.3 * min(minx, min(miny,minz));
	
	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	%[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
	[fx, fy, fz] = feval(obj.WENORK3_upwind_normal, ...
			fx, fy, fz, obj.F, fxpr, fxpl, fypf, fypb, fzpu, fzpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	fgradient = max(sqrt(fx.^2+fy.^2+fz.^2), 1e-14);
	
	nx = fx ./ fgradient;
	ny = fy ./ fgradient;
	nz = fz ./ fgradient;
	
	Sign = sign(obj.F);
	vx = Sign .* nx;
	vy = Sign .* ny;
	vz = Sign .* nz;
	
	% interpolate values for crossing points
	cpr = zeros(obj.GD3.Size, 'gpuArray');
	cpl = zeros(obj.GD3.Size, 'gpuArray');
	cpf = zeros(obj.GD3.Size, 'gpuArray');
	cpb = zeros(obj.GD3.Size, 'gpuArray');
	cpu = zeros(obj.GD3.Size, 'gpuArray');
	cpd = zeros(obj.GD3.Size, 'gpuArray');
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	step = zeros(obj.GD3.Size, 'gpuArray');
	NewA = A;
	for i=1:iter1
		NewA = ENORK2Extend(obj,NewA,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

		step = feval(obj.WENORK3_reinitialization_step, step, NewA, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Atmp1 = NewA - step;
		Atmp1 = ENORK2Extend(obj,Atmp1,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);
	
		step = feval(obj.WENORK3_reinitialization_step, step, Atmp1, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
		    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Atmp2 = Atmp1 - step;

		Atmp0_5 = 3./4. * NewA + 1./4. * Atmp2;
		Atmp0_5 = ENORK2Extend(obj,Atmp0_5,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

		step = feval(obj.WENORK3_reinitialization_step, step, Atmp0_5, mask, deltat, ...
				xpr, xpl, ypf, ypb, zpu, zpd, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
		Atmp1_5 = Atmp0_5 - step;

		NewA = 1./3. * NewA + 2./3. * Atmp1_5;
	end
	NewA = ENORK2Extend(obj,NewA,iter2,fxpr,fxpl,fypf,fypb,fzpu,fzpd,vx,vy,vz,fdeltat,cpr,cpl,cpf,cpb,cpu,cpd);

end

function NewC = ENORK2Extend(obj, C, iteration, xpr,xpl,ypf,ypb,zpu,zpd,vx,vy,vz,deltat,...
					cpr,cpl,cpf,cpb,cpu,cpd)

	NewC = C;

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.WENORK3_boundary_interpolate, ...
			cpr,cpl,cpf,cpb,cpu,cpd,xpr,xpl,ypf,ypb,zpu,zpd,NewC, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.WENORK3_extend_step,step,deltat,NewC,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		Ctmp1 = NewC - step;

		step = feval(obj.WENORK3_extend_step,step,deltat,Ctmp1,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;

		step = feval(obj.WENORK3_extend_step,step,deltat,Ctmp0_5,vx,vy,vz,...
				xpr,xpl,ypf,ypb,zpu,zpd,cpr,cpl,cpf,cpb,cpu,cpd,...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	

		Ctmp1_5 = Ctmp0_5 - step;

		%NewC = (NewC + Ctmp - step) / 2;
		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
	end
	
end
