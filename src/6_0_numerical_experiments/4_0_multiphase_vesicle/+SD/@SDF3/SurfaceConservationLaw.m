% solve the surface conservation law as a spatial conservation law for iter steps
function NewC = SurfaceConservationLaw(obj,C,vx,vy,vz,iteration,dt)
	
	NewC = C;

	step = zeros(obj.GD3.Size,'gpuArray');

	for ii=1:iteration
		step = feval(obj.spatial_finite_volume_step, step, vx, vy, vz, NewC, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp1 = NewC - step;
		Ctmp1 = obj.WENORK3Extend(Ctmp1, 10);

		step = feval(obj.spatial_finite_volume_step, step, vx, vy, vz, Ctmp1, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp2 = Ctmp1 - step;

		Ctmp0_5 = 3./4. * NewC + 1./4. * Ctmp2;
		Ctmp0_5 = obj.WENORK3Extend(Ctmp0_5, 10);

		step = feval(obj.spatial_finite_volume_step, step, vx, vy, vz, Ctmp0_5, dt, ...
				obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
				obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);
		Ctmp1_5 = Ctmp0_5 - step;

		NewC = 1./3. * NewC + 2./3. * Ctmp1_5;
		NewC = obj.WENORK3Extend(NewC, 10);
	end
end
