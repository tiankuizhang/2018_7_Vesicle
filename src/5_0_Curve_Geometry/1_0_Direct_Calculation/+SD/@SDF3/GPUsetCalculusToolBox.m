% gpu kernel version of setCalculusToolBox
function GPUsetCalculusToolBox(obj)

	obj.HPrimal = max(obj.F,0);

	obj.Fx = zeros(obj.GD3.Size,'gpuArray');
	obj.Fy = zeros(obj.GD3.Size,'gpuArray');
	obj.Fz = zeros(obj.GD3.Size,'gpuArray');
	obj.FGradMag = zeros(obj.GD3.Size,'gpuArray');

	obj.Nx = zeros(obj.GD3.Size,'gpuArray');
	obj.Ny = zeros(obj.GD3.Size,'gpuArray');
	obj.Nz = zeros(obj.GD3.Size,'gpuArray');


	obj.Fxx = zeros(obj.GD3.Size,'gpuArray');
	obj.Fyy = zeros(obj.GD3.Size,'gpuArray');
	obj.Fzz = zeros(obj.GD3.Size,'gpuArray');
	obj.Fxy = zeros(obj.GD3.Size,'gpuArray');
	obj.Fyz = zeros(obj.GD3.Size,'gpuArray');
	obj.Fzx = zeros(obj.GD3.Size,'gpuArray');
	obj.FLaplacian = zeros(obj.GD3.Size,'gpuArray');
	
	obj.MeanCurvature = zeros(obj.GD3.Size,'gpuArray');
	obj.GaussianCurvature = zeros(obj.GD3.Size,'gpuArray');

	obj.Heaviside = zeros(obj.GD3.Size,'gpuArray');
	obj.DiracDelta = zeros(obj.GD3.Size,'gpuArray');

	[obj.Fx,obj.Fy,obj.Fz,obj.FGradMag,obj.Nx,obj.Ny,obj.Nz,...
	 obj.Fxx,obj.Fyy,obj.Fzz,obj.Fxy,obj.Fyz,obj.Fzx,obj.FLaplacian,...
	 obj.MeanCurvature,obj.GaussianCurvature,...
	 obj.Heaviside,obj.DiracDelta] = feval(obj.set_calculus_toolbox, ...
		obj.Fx,obj.Fy,obj.Fz,obj.FGradMag,obj.Nx,obj.Ny,obj.Nz,...
	 	obj.Fxx,obj.Fyy,obj.Fzz,obj.Fxy,obj.Fxy,obj.Fzx,obj.FLaplacian,...
	 	obj.MeanCurvature,obj.GaussianCurvature,...
		obj.Heaviside,obj.DiracDelta, ...
		obj.F, obj.HPrimal, ...
		obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
		obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.Ds, obj.GD3.NumElt);	

end
