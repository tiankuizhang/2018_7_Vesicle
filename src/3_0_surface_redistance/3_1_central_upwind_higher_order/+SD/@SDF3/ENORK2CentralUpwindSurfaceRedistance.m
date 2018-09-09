% reinitialization with ENO and RK2 scheme for iteration times
% the boundary terms is treated the way as in the reinitialization scheme
% see 2005_Bryson_Semi-discrete central-upwind schemes with reduced dissipation for Hamilton-Jacobi equatiosn for details of contructing the numerical Hamiltonian
function NewAF = ENORK2CentralUpwindSurfaceRedistance(obj,AF,iteration)

	% this step can improve numerical accuracy 
	% note that this step can change sign of nodes and to calculate sign
	% NewAF instead of AF should be used
	% NewAF = obj.ENORK2Extend(AF,10);

	NewAF = AF;

	xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
	ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
	zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
	
	% NewAF instead of Fgpu should be used!!!!
	[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
			xpr, xpl, ypf, ypb, zpu, zpd, NewAF, obj.GD3.NumElt, ...
		    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

	% calculate the normal to the surface upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	[fx, fy, fz] = feval(obj.ENORK2_upwind_normal, ...
			fx, fy, fz, obj.F, xpr, xpl, ypf, ypb, zpu, zpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
	
	fgradient = max( sqrt(fx.^2+fy.^2+fz.^2), 1e-14);
	
	nx = fx ./ fgradient;
	ny = fy ./ fgradient;
	nz = fz ./ fgradient;

	% calculate the sign of the original auxiliary level set function
	% NewAF instead of AF should be used
	Sign = zeros(obj.GD3.Size, 'gpuArray');
	Sign(NewAF>0) = 1.;
	Sign(NewAF<0) = -1.;
	
	% CFL condition imposed on time step
	minx = min(xpr,xpl);
	miny = min(ypf,ypb);
	minz = min(zpu,zpd);
	deltat = 0.3 * min(minx, min(miny,minz));

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.ENORK2_surface_redistance_step,step,NewAF,Sign,deltat,nx,ny,nz,...
			   		xpr,xpl,ypf,ypb,zpu,zpd,...	
			    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
					obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		AFtmp = NewAF - step;
		step = feval(obj.ENORK2_surface_redistance_step,step,AFtmp,Sign,deltat,nx,ny,nz,...
			   		xpr,xpl,ypf,ypb,zpu,zpd,...	
			    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
					obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
		NewAF = (NewAF + AFtmp - step) / 2;
	end

	%NewAF = obj.ENORK2Extend(NewAF,50);

end
