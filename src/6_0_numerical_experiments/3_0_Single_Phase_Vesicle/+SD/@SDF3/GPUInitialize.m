% initialize GPU functions
function GPUInitialize(obj)

	%fprintf('\t\t creating CUDAKernel ...')

	%obj.ThreadBlockSize = [obj.GD3.mrows,floor(512/obj.GD3.mrows),1];
	obj.ThreadBlockSize = [obj.GD3.mrows,floor(256/obj.GD3.mrows),1];
	obj.GridSize = [ceil(obj.GD3.mrows/obj.ThreadBlockSize(1)), ...
					ceil(obj.GD3.ncols/obj.ThreadBlockSize(2)), ...
					ceil(obj.GD3.lshts/obj.ThreadBlockSize(3))];

	% functions used by 3rd order reinitialization scheme and other schemes
	[status, cmdout] = system('make -C CUDA_Code/1_0_ENORK2_Reinitialization');
	if status, error(cmdout), end % if compilation failed, throw and output error message

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

	% function used by 4th order WENO reinitialization scheme
	[status, cmdout] = system('make -C CUDA_Code/1_1_WENORK3_Reinitialization');
	if status, error(cmdout), end % if compilation failed, throw and output error message

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

	% function used by 6th order WENO5 reinitialization scheme
	[status, cmdout] = system('make -C CUDA_Code/1_3_WENO5RK3_Reinitialization');
	if status, error(cmdout), end % if compilation failed, throw and output error message

	obj.sixth_boundary_correction = parallel.gpu.CUDAKernel('CUDA_Code/1_3_WENO5RK3_Reinitialization/boundary_location.ptx', ...
															'CUDA_Code/1_3_WENO5RK3_Reinitialization/boundary_location.cu', ...
															'boundary_location');
	obj.sixth_boundary_correction.ThreadBlockSize = obj.ThreadBlockSize;
	obj.sixth_boundary_correction.GridSize = obj.GridSize;

	obj.WENO5RK3_reinitialization_step = parallel.gpu.CUDAKernel(	'CUDA_Code/1_3_WENO5RK3_Reinitialization/re_step.ptx', ...
																	'CUDA_Code/1_3_WENO5RK3_Reinitialization/re_step.cu', ...
																	're_step');
	obj.WENO5RK3_reinitialization_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENO5RK3_reinitialization_step.GridSize = obj.GridSize;

	% functions used by the 3rd order extend scheme and other schemes
	[status, cmdout] = system('make -C CUDA_Code/2_0_ENORK2_Extend');
	if status, error(cmdout), end % if compilation failed, throw and output error message

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

	% functions used by the 4th order WENO extend scheme and other schemes
	[status, cmdout] = system('make -C CUDA_Code/2_1_WENORK3_Extend');
	if status, error(cmdout), end % if compilation failed, throw and output error message

	obj.WENORK3_upwind_normal = parallel.gpu.CUDAKernel('CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
														'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
														'upwind_normal');
	obj.WENORK3_upwind_normal.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_upwind_normal.GridSize = obj.GridSize;

	obj.WENORK3_extend_step = parallel.gpu.CUDAKernel(  'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
														'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
														'extend_step');
	obj.WENORK3_extend_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_extend_step.GridSize = obj.GridSize;

	obj.WENORK3_boundary_interpolate = parallel.gpu.CUDAKernel( 'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.ptx', ...
																'CUDA_Code/2_1_WENORK3_Extend/wenork3_extend.cu', ...
																'boundary_interpolate');
	obj.WENORK3_boundary_interpolate.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENORK3_boundary_interpolate.GridSize = obj.GridSize;		

	% functions used by 6th order WENO extend scheme and other schemes
	[status, cmdout] = system('make -C CUDA_Code/2_2_WENO5RK3_Extend');
	if status, error(cmdout), end;

	obj.WENO5RK3_upwind_normal = parallel.gpu.CUDAKernel('CUDA_Code/2_2_WENO5RK3_Extend/upwind_normal.ptx', ...
														 'CUDA_Code/2_2_WENO5RK3_Extend/upwind_normal.cu', ...
														 'upwind_normal');
	obj.WENO5RK3_upwind_normal.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENO5RK3_upwind_normal.GridSize = obj.GridSize;

	obj.WENO5RK3_extend_step = parallel.gpu.CUDAKernel( 'CUDA_Code/2_2_WENO5RK3_Extend/extend_step.ptx', ...
														'CUDA_Code/2_2_WENO5RK3_Extend/extend_step.cu', ...
														'extend_step');
	obj.WENO5RK3_extend_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENO5RK3_extend_step.GridSize = obj.GridSize;

	obj.WENO5RK3_boundary_interpolant = parallel.gpu.CUDAKernel('CUDA_Code/2_2_WENO5RK3_Extend/boundary_interpolant.ptx', ...
																'CUDA_Code/2_2_WENO5RK3_Extend/boundary_interpolant.cu', ...
																'boundary_interpolant');
	obj.WENO5RK3_boundary_interpolant.ThreadBlockSize = obj.ThreadBlockSize;
	obj.WENO5RK3_boundary_interpolant.GridSize = obj.GridSize;

	% functions used by GPUsetCalculusTooBox scheme
	[status, cmdout] = system('make -C CUDA_Code/4_0_Calculus_ToolBox');
	if status, error(cmdout), end % if compilation failed, throw and output error message

	obj.set_calculus_toolbox = parallel.gpu.CUDAKernel( 'CUDA_Code/4_0_Calculus_ToolBox/set_calculus_toolbox.ptx', ...
														'CUDA_Code/4_0_Calculus_ToolBox/set_calculus_toolbox.cu', ...
														'set_calculus_toolbox');
	obj.set_calculus_toolbox.ThreadBlockSize = obj.ThreadBlockSize;
	obj.set_calculus_toolbox.GridSize = obj.GridSize;

	obj.auxi_set_calculus_toolbox = parallel.gpu.CUDAKernel('CUDA_Code/4_0_Calculus_ToolBox/auxi_set_calculus_toolbox.ptx', ...
															'CUDA_Code/4_0_Calculus_ToolBox/auxi_set_calculus_toolbox.cu', ...
															'auxi_set_calculus_toolbox');
	obj.auxi_set_calculus_toolbox.ThreadBlockSize = obj.ThreadBlockSize;
	obj.auxi_set_calculus_toolbox.GridSize = obj.GridSize;

	obj.upwind_derivative = parallel.gpu.CUDAKernel('CUDA_Code/4_0_Calculus_ToolBox/upwind_derivative.ptx', ...
													'CUDA_Code/4_0_Calculus_ToolBox/upwind_derivative.cu', ...
													'upwind_derivative');
	obj.upwind_derivative.ThreadBlockSize = obj.ThreadBlockSize;
	obj.upwind_derivative.GridSize = obj.GridSize;

	obj.ga_set_calculus_toolbox = parallel.gpu.CUDAKernel(  'CUDA_Code/4_0_Calculus_ToolBox/ga_set_calculus_toolbox.ptx', ...
															'CUDA_Code/4_0_Calculus_ToolBox/ga_set_calculus_toolbox.cu', ...
															'ga_set_calculus_toolbox');
	obj.ga_set_calculus_toolbox.ThreadBlockSize = obj.ThreadBlockSize;
	obj.ga_set_calculus_toolbox.GridSize = obj.GridSize;

	obj.advection_step = parallel.gpu.CUDAKernel('CUDA_Code/4_0_Calculus_ToolBox/advection_step.ptx', ...
												 'CUDA_Code/4_0_Calculus_ToolBox/advection_step.cu', ...
												 'advection_step');
	obj.advection_step.ThreadBlockSize = obj.ThreadBlockSize;
	obj.advection_step.GridSize = obj.GridSize;

	% functions used for calculating weno derivatives
	[status, cmdout] = system('make -C CUDA_Code/5_0_Nonoscillating_Reconstruction');
	if status, error(cmdout), end % if compilation failed, throw and output error message

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

	%fprintf('finished !\n')
end








