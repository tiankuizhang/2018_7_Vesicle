all: upwind_normal.ptx boundary_interpolant.ptx extend_step.ptx

upwind_normal.ptx boundary_interpolant.ptx extend_step.ptx: shared_utilities.cup shared_utilities.cuh

upwind_normal.ptx: upwind_normal.cu 
	nvcc -ptx upwind_normal.cu -o upwind_normal.ptx

boundary_interpolant.ptx: boundary_interpolant.cu 
	nvcc -ptx boundary_interpolant.cu -o boundary_interpolant.ptx

extend_step.ptx: extend_step.cu 
	nvcc -ptx extend_step.cu -o extend_step.ptx
