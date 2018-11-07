% initialize GPU functions
function GPUInitialize(obj)
	
	fprintf("start compiling .cu into .ptx files . . .")

	obj.ThreadBlockSize = [obj.GD3.mrows,floor(512/obj.GD3.mrows),1];
	obj.GridSize = [ceil(obj.GD3.mrows/obj.ThreadBlockSize(1)), ...
					ceil(obj.GD3.ncols/obj.ThreadBlockSize(2)), ...
					ceil(obj.GD3.lshts/obj.ThreadBlockSize(3))];

	% functions used by reinitialization scheme and other schemes
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu -o CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx');

	obj.ENORK2_boundary_correction = parallel.gpu.CUDAKernel('CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx', ...
															 'CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu', ...
															 'boundary_correction');
	obj.ENORK2_boundary_correction.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ENORK2_boundary_correction.GridSize = obj.GridSize;

	obj.ENORK2_reinitiliaztion_step = parallel.gpu.CUDAKernel('CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx', ...
															  'CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu', ...
															  're_step');
	obj.ENORK2_reinitiliaztion_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ENORK2_reinitiliaztion_step.GridSize = obj.GridSize;

	% function used by WENO reinitialization scheme
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.cu -o CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.ptx');

	obj.cubic_boundary_correction = parallel.gpu.CUDAKernel('CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.ptx', ...
															'CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.cu', ...
															'boundary_correction');
	obj.cubic_boundary_correction.ThreadBlockSize = obj.ThreadBlockSize;
	obj.cubic_boundary_correction.GridSize = obj.GridSize;

	obj.WENORK3_reinitialization_step = parallel.gpu.CUDAKernel('CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.ptx', ...
																'CUDA_Code/1_1_WENORK3_Reinitialization/weno_rk3_reinitialization.cu', ...
																're_step');
	obj.WENORK3_reinitialization_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_reinitialization_step.GridSize = obj.GridSize;

	% function used by WENO5 reinitialization scheme
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.cu -o CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.ptx');

	obj.sixth_boundary_correction = parallel.gpu.CUDAKernel('CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.ptx', ...
															'CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.cu', ...
															'boundary_correction');
	obj.sixth_boundary_correction.ThreadBlockSize = obj.ThreadBlockSize;
	obj.sixth_boundary_correction.GridSize = obj.GridSize;

	obj.WENO5RK3_reinitialization_step = parallel.gpu.CUDAKernel(	'CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.ptx', ...
																	'CUDA_Code/1_3_WENO5RK3_Reinitialization/weno5_rk3_reinitialization.cu', ...
																	're_step');
	obj.WENO5RK3_reinitialization_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENO5RK3_reinitialization_step.GridSize = obj.GridSize;

	% functions used by the extend scheme and other schemes
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu -o CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx');

	obj.ENORK2_upwind_normal = parallel.gpu.CUDAKernel(	'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx', ...
		   												'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu', ...
														'upwind_normal');
	obj.ENORK2_upwind_normal.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ENORK2_upwind_normal.GridSize = obj.GridSize;		

	obj.ENORK2_extend_step = parallel.gpu.CUDAKernel('CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx', ...
		   											 'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu', ...
												   	 'extend_step');
	obj.ENORK2_extend_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ENORK2_extend_step.GridSize = obj.GridSize;		

	obj.ENORK2_boundary_interpolate = parallel.gpu.CUDAKernel('CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx', ...
		   											 		  'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu', ...
		          										   	  'boundary_interpolate');
	obj.ENORK2_boundary_interpolate.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ENORK2_boundary_interpolate.GridSize = obj.GridSize;		

	% functions used by the WENO extend scheme and other schemes
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu -o CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx');

	obj.WENORK3_upwind_normal = parallel.gpu.CUDAKernel('CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
		   											    'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
												   	    'upwind_normal');
	obj.WENORK3_upwind_normal.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_upwind_normal.GridSize = obj.GridSize;		

	obj.WENORK3_extend_step = parallel.gpu.CUDAKernel('CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
		   											  'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
												   	  'extend_step');
	obj.WENORK3_extend_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_extend_step.GridSize = obj.GridSize;		

	obj.WENORK3_boundary_interpolate = parallel.gpu.CUDAKernel('CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
		   											 		   'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
		          										   	   'boundary_interpolate');
	obj.WENORK3_boundary_interpolate.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_boundary_interpolate.GridSize = obj.GridSize;		

	% functions used by GPUsetCalculusTooBox scheme
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.cu -o CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.ptx');

	obj.set_calculus_toolbox = parallel.gpu.CUDAKernel('CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.ptx', ...
													   'CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.cu', ...
													   'set_calculus_toolbox');
	obj.set_calculus_toolbox.ThreadBlockSize = obj.ThreadBlockSize;
	obj.set_calculus_toolbox.GridSize = obj.GridSize;

	obj.auxi_set_calculus_toolbox = parallel.gpu.CUDAKernel('CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.ptx', ...
													   'CUDA_Code/4_0_Calculus_ToolBox/calculus_toolbox.cu', ...
													   'auxi_set_calculus_toolbox');
	obj.auxi_set_calculus_toolbox.ThreadBlockSize = obj.ThreadBlockSize;
	obj.auxi_set_calculus_toolbox.GridSize = obj.GridSize;

	% functions used for calculating weno derivatives
	system('nvcc -Wno-deprecated-gpu-targets -ptx CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.cu -o CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.ptx');

	obj.GD3.weno_derivative = parallel.gpu.CUDAKernel('CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.ptx', ...
												  	  'CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.cu', ...
												  	  'weno_derivative');
	obj.GD3.weno_derivative.ThreadBlockSize = obj.ThreadBlockSize;
	obj.GD3.weno_derivative.GridSize = obj.GridSize;

	obj.surface_conservation_step = parallel.gpu.CUDAKernel('CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.ptx', ...
												  			'CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.cu', ...
												  	  		'surface_conservation_step');
	obj.surface_conservation_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.surface_conservation_step.GridSize = obj.GridSize;

	obj.spatial_finite_volume_step = parallel.gpu.CUDAKernel('CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.ptx', ...
												  			'CUDA_Code/5_0_Nonoscillating_Reconstruction/nonoscillating_interpolant.cu', ...
												  	  		'spatial_finite_volume_step');
	obj.spatial_finite_volume_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.spatial_finite_volume_step.GridSize = obj.GridSize;

	fprintf("finished compilation!\n")
end








