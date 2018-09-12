% calculate WENO derivatives of Field with CUDA kernel function
function [WENO_back_x, WENO_fore_x, WENO_back_y, WENO_fore_y, WENO_back_z, WENO_fore_z] = GPUWENODerivative(obj, Field)

	WENO_back_x = zeros(obj.Size, 'gpuArray');
	WENO_fore_x = zeros(obj.Size, 'gpuArray');
	WENO_back_y = zeros(obj.Size, 'gpuArray');
	WENO_fore_y = zeros(obj.Size, 'gpuArray');
	WENO_back_z = zeros(obj.Size, 'gpuArray');
	WENO_fore_z = zeros(obj.Size, 'gpuArray');

	[WENO_back_x, WENO_fore_x, WENO_back_y, WENO_fore_y, WENO_back_z, WENO_fore_z] = ...
		  feval(obj.weno_derivative, ...
				WENO_back_x, WENO_fore_x, WENO_back_y, WENO_fore_y, WENO_back_z, WENO_fore_z,...
			   	Field, obj.mrows, obj.ncols, obj.lshts, obj.Dx, obj.Dy, obj.Dz);	
end
