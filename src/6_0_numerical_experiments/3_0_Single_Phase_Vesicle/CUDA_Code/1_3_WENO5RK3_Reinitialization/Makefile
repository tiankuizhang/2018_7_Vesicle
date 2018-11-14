all: re_step.ptx boundary_location.ptx

re_step.ptx: re_step.cu shared_utilities.hpp
	nvcc -ptx re_step.cu -o re_step.ptx

boundary_location.ptx: boundary_location.cu shared_utilities.hpp
	nvcc -ptx boundary_location.cu -o boundary_location.ptx
