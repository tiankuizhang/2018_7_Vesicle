/*******************************************************************************
 * given a velocity field (vx,vy,vz) and a scalar field lsf
 * calculate v*grad(lsf) with WENO scheme
 ******************************************************************************/
#include "shared_utilities.cuh"
#include "shared_utilities.cup"

// calculate forward and back weno derivative in the direction specified by step[i]
// i is 0 for x, 1 for y, 2 for z
__device__ inline
void weno_derivative(double & df, double & db, double const * lsf, int const step[3], int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges, double ds)
{
	int stencil[7]; // index of the 7 point stencil
	for(int i=0; i<7; i++){
		int row_shift = row_idx + (i - 3) * step[0];
		int col_shift = col_idx + (i - 3) * step[1];
		int pge_shift = pge_idx + (i - 3) * step[2];
		stencil[i] = sub2ind(row_shift, col_shift, pge_shift, rows, cols, pges);
	}
	
	double v[6]; // cell average of derivatives
	for(int i=0; i<6; i++){
		v[i] = (lsf[ stencil[i+1] ] - lsf[ stencil[i] ]) / ds;
	}

	db = weno_onesided_derivative(v[0],v[1],v[2],v[3],v[4]);
	df = weno_onesided_derivative(v[5],v[4],v[3],v[2],v[1]);
}

__global__
void advection_step(double * astep, double const * vx, double const * vy, double const * vz, double const * lsf, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int compass[3][3] = {{0,1,0},{1,0,0},{0,0,1}};

	double xL,xR,yB,yF,zD,zU;

	weno_derivative(xR,xL,lsf,compass[0],row_idx,col_idx,pge_idx,rows,cols,pges,dx);
	weno_derivative(yF,yB,lsf,compass[1],row_idx,col_idx,pge_idx,rows,cols,pges,dy);
	weno_derivative(zU,zD,lsf,compass[2],row_idx,col_idx,pge_idx,rows,cols,pges,dz);

	astep[ind] = (min2(0.0,vx[ind]) * xR + max2(0.0,vx[ind]) * xL +
				  min2(0.0,vy[ind]) * yF + max2(0.0,vy[ind]) * yB +
				  min2(0.0,vz[ind]) * zU + max2(0.0,vz[ind]) * zD );


}




































