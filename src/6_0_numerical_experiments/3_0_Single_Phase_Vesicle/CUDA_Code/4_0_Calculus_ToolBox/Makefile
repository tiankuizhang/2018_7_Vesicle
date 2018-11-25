all: set_calculus_toolbox.ptx auxi_set_calculus_toolbox.ptx upwind_derivative.ptx

set_calculus_toolbox.ptx auxi_set_calculus_toolbox.ptx upwind_derivative.ptx: shared_utilities.cup shared_utilities.cuh 

set_calculus_toolbox.ptx: set_calculus_toolbox.cu
	nvcc -ptx set_calculus_toolbox.cu -o set_calculus_toolbox.ptx

auxi_set_calculus_toolbox.ptx: auxi_set_calculus_toolbox.cu
	nvcc -ptx auxi_set_calculus_toolbox.cu -o auxi_set_calculus_toolbox.ptx

upwind_derivative.ptx: upwind_derivative.cu
	nvcc -ptx upwind_derivative.cu -o upwind_derivative.ptx
