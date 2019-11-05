% solve surface conservation law for iter steps
function NewC = WENORK3SurfaceConservationLaw(obj,C,vx,vy,vz,iteration,dt)

	NewC = C;

	vd = obj.SurfaceDivergence(vx,vy,vz);

	step = zeros(obj.GD3.Size, 'gpuArray');

	for i=1:iteration
		step = feval(obj.surface_conservation_step, step, vx, vy, vz, vd, NewC, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp1 = NewC - step;

		step = feval(obj.surface_conservation_step, step, vx, vy, vz, vd, Ctmp1, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;

		step = feval(obj.surface_conservation_step, step, vx, vy, vz, vd, Ctmp0_5, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp1_5 = Ctmp0_5 - step;

		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
	end


end
