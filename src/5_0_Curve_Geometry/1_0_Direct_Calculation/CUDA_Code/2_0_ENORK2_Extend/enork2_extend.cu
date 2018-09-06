/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the extension scheme 
 ******************************************************************************/

typedef struct
{
	double sR;
	double sL;
} double_eno_derivative;

__device__ inline
double max2(double x, double y)
{
	return (x<y) ? y : x;
}

__device__ inline
double min2(double x, double y)
{
	return (x<y) ? x : y;
}

__device__ inline
double min_mod(double x, double y)
{
	return (x*y<0) ? 0.0 : (fabs(x)<fabs(y) ? x : y);
}


__device__ inline
double sign(double x)
{
	return (x>0) ? 1.0 : -1.0;
}


// convert subindex to linear index
// periodic boundary conditions are assumed
__device__ inline
int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges)
{	
	int row_idxn = min2(rows-1, max2(0, row_idx));
	int col_idxn = min2(cols-1, max2(0, col_idx));
	int pge_idxn = min2(pges-1, max2(0, pge_idx));

	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

/****************************************************************************** 
 * calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
 * 
 ******************************************************************************/
__device__ inline
double_eno_derivative eno_derivative_field( double v4, double v1, double v0, double v2, double v3, double ds, double dis_f, double dis_b, double v_forward, double v_backward)
{
	double p2m;
	double_eno_derivative eno_d;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	double v_f = (dis_f==ds) ? v2 : v_forward;
	eno_d.sR = (v_f - v0) / dis_f - dis_f * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	double v_b = (dis_b==ds) ? v1 : v_backward;
	eno_d.sL = (v0 - v_b) / dis_b + dis_b * p2m;

	return eno_d;

}

/*******************************************************************************
 * calculate upwind normal with ENO scheme at a single point
 * along a single direction
 *******************************************************************************/                                                                                 
__device__ inline
double upwind_normal_point( double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	double vr = (pr==ds) ? v2 : 0;
	double sR = (vr - v0) / pr - pr * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	double vl = (pl==ds) ? v1 : 0;
	double sL = (v0 - vl) / pl + pl * p2m;

	return (fabs(vr) < fabs(vl)) ? sR : sL;
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

	int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);

	nx[ind] = upwind_normal_point( lsf[left2], lsf[left], lsf[ind], lsf[right], lsf[right2], xpr[ind], xpl[ind], dx);

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	ny[ind] = upwind_normal_point( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	nz[ind] = upwind_normal_point( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);

}	

// modify forward/backward(c_f/b) values at v0:[v4,v1,v0,v2,v3]
// if dis_b/f!=ds, then there is boundary nearby, backward/forward ENO is then constructed
// C(x) = c2*x^2 + c1*x + c0, c(-ds/2) = v0, c(ds/2) = v2 assuming boundary is between v0,v2
// and used to calculate c_b/f at boundary crossing nodes
__device__ inline
void cubic_eno_interp(double & c_forward, double & c_backward, double dis_f, double dis_b, double ds, double v4, double v1, double v0, double v2, double v3)
{
	// if there is a boundary in the forward direction
	c_forward = 0;
	c_backward = 0;
	if(dis_f!=ds){
		double p2  = v2 - 2.0 * v0 + v1;
		double p2r = v3 - 2.0 * v2 + v0;
		double c2 = 0.5 * min_mod(p2,p2r) / pow(ds,2); // coefficient for second order term
		double c1 = (v2 - v0) / ds; // coefficient for the linear term
		double c0 = (v2 + v0) / 2 - c2 * pow(ds,2) / 4; // constant term
		double x = dis_f - ds / 2; // coordinate for the boundary point
		c_forward = c2 * x * x + c1 * x + c0; // interpolated value at the right boundary
		// note that the formula works even if c2 is very small
	    // then we have a linear interpolation
	}
	// if there is a boundary in the backward direction
	if(dis_b!=ds){
		double p2  = v2 - 2.0 * v0 + v1;
		double p2l = v0 - 2.0 * v1 + v4;
		double c2 = 0.5 * min_mod(p2,p2l) / pow(ds,2);
		double c1 = (v0 - v1) / ds;
		double c0 = (v1 + v0) / 2 - c2 * pow(ds,2) / 4;
		double x = ds / 2 - dis_b;
		c_backward = c2 * x * x + c1 * x + c0;
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

	int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);

	cubic_eno_interp(cpr[ind],cpl[ind],xpr[ind],xpl[ind],dx,lsf[left2],lsf[left],lsf[ind],lsf[right],lsf[right2]);

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	cubic_eno_interp(cpf[ind],cpb[ind],ypf[ind],ypb[ind],dy,lsf[back2],lsf[back],lsf[ind],lsf[front],lsf[front2]);

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	cubic_eno_interp(cpu[ind],cpd[ind],zpu[ind],zpd[ind],dz,lsf[down2],lsf[down],lsf[ind],lsf[up],lsf[up2]);

}


// calculate extend step
// now lsf represents a scalar field (not the level set function)
__global__
void extend_step(double * step, double const * deltat, double const * lsf, double const * vx, double const * vy, double const * vz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, double const * cpr, double const * cpl, double const * cpf, double const * cpb, double const * cpu, double const * cpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dx = eno_derivative_field( lsf[left2], lsf[left], lsf[ind], lsf[right], lsf[right2], dx, xpr[ind], xpl[ind], cpr[ind], cpl[ind]);
	double xR = eno_dx.sR;
	double xL = eno_dx.sL;


	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative_field( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], dy, ypf[ind], ypb[ind], cpf[ind], cpb[ind]);
	double yF = eno_dy.sR;
	double yB = eno_dy.sL;


	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative_field( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], dz, zpu[ind], zpd[ind], cpu[ind], cpd[ind]);
	double zU = eno_dz.sR;
	double zD = eno_dz.sL;

	step[ind] = (min2(0,vx[ind]) * xR + max2(0,vx[ind]) * xL + 
				 min2(0,vy[ind]) * yF + max2(0,vy[ind]) * yB + 
				 min2(0,vz[ind]) * zU + max2(0,vz[ind]) * zD ) * deltat[ind];

}








