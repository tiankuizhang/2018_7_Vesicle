/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the extension scheme 
 ******************************************************************************/
#include "shared_utilities.cuh"
#include "shared_utilities.cup"

// suppose that (3*s,v3) is a boundary node and we know the distance from this node to the boundary
// sixth_interp will interpolate values onto the boundary node from the stencil
// (0,v0),(s,v1),(2s,v2),(3s,v3),(4s,v4),(5s,v5) if p2*p3<0 for c_backward
// (s,v1),(2s,v2),(3s,v3),(4s,v4),(5s,v5),(6s,v6) if p3*p4<0 for c_forward
__device__ inline
void sixth_interp(double & c_forward, double & c_backward, double dis_f, double dis_b, double s, double v0, double v1, double v2, double v3, double v4, double v5, double v6)
{
	c_forward = 0.0;
	c_backward = 0.0;

	// if there is a boundary in the forward direction
	if(dis_f<s){
		double xf[6] = {s,2*s,3*s,4*s,5*s,6*s};
		double yf[6] = {v1,v2,v3,v4,v5,v6};
		Poly<5> p5_f;
		p5_f.setInterpolation(xf,yf);
		c_forward = p5_f(3*s + dis_f);
	}
	// if there is a boundary in the backward direction
	if(dis_f<s){
		double xf[6] = {0,s,2*s,3*s,4*s,5*s};
		double yf[6] = {v0,v1,v2,v3,v4,v5};
		Poly<5> p5_b;
		p5_b.setInterpolation(xf,yf);
		c_backward = p5_b(3*s - dis_b);
	}
}


// interpolate values at boundary points
__global__
void boundary_interpolate(double * cpr, double * cpl, double * cpf, double * cpb, double * cpu, double * cpd, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, double const * lsf, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int right  = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int right3 = sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left   = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2  = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3  = sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

	sixth_interp(cpr[ind],cpl[ind],xpr[ind],xpl[ind],dx,lsf[left3],lsf[left2],lsf[left],lsf[ind],lsf[right],lsf[right2],lsf[right3]);

	int front  = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int front3 = sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);
	int back   = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2  = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3  = sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

	sixth_interp(cpf[ind],cpb[ind],ypf[ind],ypb[ind],dy,lsf[back3],lsf[back2],lsf[back],lsf[ind],lsf[front],lsf[front2],lsf[front3]);

	int up    = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	int up2   = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
	int up3   = sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);
	int down  = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 = sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

	sixth_interp(cpu[ind],cpd[ind],zpu[ind],zpd[ind],dz,lsf[down3],lsf[down2],lsf[down],lsf[ind],lsf[up],lsf[up2],lsf[up3]);
}







