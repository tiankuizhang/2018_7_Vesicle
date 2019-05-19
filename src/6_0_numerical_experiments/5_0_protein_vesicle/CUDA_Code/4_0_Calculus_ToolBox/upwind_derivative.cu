// uwpind_derivative will calculate gradient of a field living on the
// 2d surface embedded in the 3d space
// due to the embedding nature, we may use geometrical information in the
// level set to avoid interpolation through discontinuities
// the gradient will be the gradient selected in a normal flow from the boundary

#include "shared_utilities.cuh"
#include "shared_utilities.cup"

__global__
void upwind_derivative(double * fx, double * fy, double * fz, double const * lsf, double const * vx, double const * vy, double const * vz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, double const * cpr, double const * cpl, double const * cpf, double const * cpb, double const * cpu, double const * cpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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
	double v_fore, v_back;

	p4 = lsf[ind];

	int rght1 = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int rght2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int rght3 = sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left1 = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3 = sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

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

	v_fore = cpr[ind];
	v_back = cpl[ind];

	double xR, xL;
	weno_derivative_boundary(xR,xL,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dx,v_fore,v_back);

	fx[ind] = (vx[ind]>0) ? xL : ( (vx[ind]<0) ? xR : (xL+xR)/2.0 );

	int frnt1 = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	int frnt2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int frnt3 = sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);
	int back1 = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3 = sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

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

	v_fore = cpf[ind];
	v_back = cpb[ind];

	double yF, yB;
	weno_derivative_boundary(yF,yB,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dy,v_fore,v_back);

	fy[ind] = (vy[ind]>0) ? yB : ( (vx[ind]<0) ? yF : (yB+yF)/2.0 );

	int upup1 = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	int upup2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
	int upup3 = sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);
	int down1 = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 = sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

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

	v_fore = cpu[ind];
	v_back = cpd[ind];

	double zU, zD;
	weno_derivative_boundary(zU,zD,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dz,v_fore,v_back);

	fz[ind] = (vz[ind]>0) ? zD : ( (vz[ind]<0) ? zU : (zD+zU)/2.0 );

}










