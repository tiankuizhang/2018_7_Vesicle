#include "shared_utilities.cuh"
#include "shared_utilities.cup"

__device__ inline
double upwind_normal_point(double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
	double d_back, d_fore;
	weno_derivative_boundary(d_fore,d_back,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds,0.0,0.0);

	return (fabs(p5)<fabs(p3)) ? d_fore : d_back;
}

// calculate the upwind normal
__global__
void upwind_normal(double * nx, double * ny, double * nz, double const * lsf, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	double p1,p2,p3,p4,p5,p6,p7;
	double r1,r2,r3,l1,l2,l3;

	p4 = lsf[ind];

	int rght1 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int rght2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int rght3 	= sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left1 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3 	= sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

	p1 = lsf[left3];
	p2 = lsf[left2];
	p3 = lsf[left1];
	p5 = lsf[rght1];
	p6 = lsf[rght2];
	p7 = lsf[rght3];

	r1 = xpr[ind];
	r2 = xpr[rght1];
	r3 = xpr[rght2];

	l1 = xpl[ind];
	l2 = xpl[left1];
	l3 = xpl[left2];

	nx[ind] = upwind_normal_point(p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dx);

	int frnt1 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int frnt2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);	
	int frnt3 	= sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);	
	int back1 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3 	= sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

	p1 = lsf[back3];
	p2 = lsf[back2];
	p3 = lsf[back1];
	p5 = lsf[frnt1];
	p6 = lsf[frnt2];
	p7 = lsf[frnt3];

	r1 = ypf[ind];
	r2 = ypf[frnt1];
	r3 = ypf[frnt2];

	l1 = ypb[ind];
	l2 = ypb[back1];
	l3 = ypb[back2];

	ny[ind] = upwind_normal_point(p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dy);

	int upup1	= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int upup2	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int upup3	= sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);	
	int down1 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 	= sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

	p1 = lsf[down3];
	p2 = lsf[down2];
	p3 = lsf[down1];
	p5 = lsf[upup1];
	p6 = lsf[upup2];
	p7 = lsf[upup3];

	r1 = zpu[ind];
	r2 = zpu[upup1];
	r3 = zpu[upup2];

	l1 = zpd[ind];
	l2 = zpd[down1];
	l3 = zpd[down2];

	nz[ind] = upwind_normal_point(p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dz);
	
}
