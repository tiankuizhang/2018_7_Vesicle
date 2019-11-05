/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to calculate 
 * weno derivatives
 ******************************************************************************/

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

// convert subindex to linear index
__device__ inline
int sub2ind(int row_idx, int col_idx, int pge_idx, int rows, int cols, int pges)
{	
	int row_idxn = min2(rows-1, max2(0, row_idx));
	int col_idxn = min2(cols-1, max2(0, col_idx));
	int pge_idxn = min2(pges-1, max2(0, pge_idx));

	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

__device__ inline
double weno_onesided_derivative(double v1, double v2, double v3, double v4, double v5)
{
	// different choices of ENO derivatives
	double phi1 = 1./3. * v1 	- 7./6. * v2 	+ 11./6. * v3;
	double phi2 = -1./6. * v2 	+ 5./6. * v3 	+ 1./3. * v4;
	double phi3 = 1./3. * v3	+ 5./6. * v4	- 1./6. * v5;

	// smoothness parameter
	double S1 = 13./12. * pow((v1 - 2*v2 + v3),2) + 1./4. * pow((v1 - 4*v2 + 3*v3),2);
	double S2 = 13./12. * pow((v2 - 2*v3 + v4),2) + 1./4. * pow((v2 - v4),2);
	double S3 = 13./12. * pow((v3 - 2*v4 + v5),2) + 1./4. * pow((3*v3 - 4*v4 + v5),2);

	double epsilon = 1e-6;
	double alpha1 = 0.1 / pow( (S1 + epsilon), 2);
	double alpha2 = 0.6 / pow( (S2 + epsilon), 2);
	double alpha3 = 0.3 / pow( (S3 + epsilon), 2);

	// weights for each stencil
	double sum = alpha1 + alpha2 + alpha3;
	double omega1 = alpha1 / sum;
	double omega2 = alpha2 / sum;
	double omega3 = alpha3 / sum;

	return (omega1*phi1 + omega2*phi2 + omega3*phi3);

}

__global__ 
void weno_derivative(double * WENO_back_x, double * WENO_fore_x, double * WENO_back_y, double * WENO_fore_y, double * WENO_back_z, double * WENO_fore_z, double const * lsf, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int rght1 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int rght2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int rght3 	= sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left1 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3 	= sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

	double v1 = (lsf[left2] - lsf[left3]) / dx;
	double v2 = (lsf[left1] - lsf[left2]) / dx;
	double v3 = (lsf[ind]   - lsf[left1]) / dx;
	double v4 = (lsf[rght1] - lsf[ind])   / dx;
	double v5 = (lsf[rght2] - lsf[rght1]) / dx;
	double v6 = (lsf[rght3] - lsf[rght2]) / dx;

	WENO_back_x[ind] = weno_onesided_derivative(v1,v2,v3,v4,v5);
	WENO_fore_x[ind] = weno_onesided_derivative(v6,v5,v4,v3,v2);

	int frnt1 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int frnt2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);	
	int frnt3 	= sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);	
	int back1 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3 	= sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

	v1 = (lsf[back2] - lsf[back3]) / dy;
	v2 = (lsf[back1] - lsf[back2]) / dy;
	v3 = (lsf[ind]	 - lsf[back1]) / dy;
	v4 = (lsf[frnt1] - lsf[ind])  / dy;
	v5 = (lsf[frnt2] - lsf[frnt1]) / dy;
	v6 = (lsf[frnt3] - lsf[frnt2]) / dy;
	
	WENO_back_y[ind] = weno_onesided_derivative(v1,v2,v3,v4,v5);
	WENO_fore_y[ind] = weno_onesided_derivative(v6,v5,v4,v3,v2);

	int upup1	= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int upup2	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int upup3	= sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);	
	int down1 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 	= sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

	v1 = (lsf[down2] - lsf[down3]) / dz;
	v2 = (lsf[down1] - lsf[down2]) / dz;
	v3 = (lsf[ind]	 - lsf[down1]) / dz;
	v4 = (lsf[upup1] - lsf[ind])   / dz;
	v5 = (lsf[upup2] - lsf[upup1]) / dz;
	v6 = (lsf[upup3] - lsf[upup2]) / dz;

	WENO_back_z[ind] = weno_onesided_derivative(v1,v2,v3,v4,v5);
	WENO_fore_z[ind] = weno_onesided_derivative(v6,v5,v4,v3,v2);
}	


// calculate numerical Hamiltonian for surface conservation law equation with input vx,vy,vz and
// surface divergence vd of v
__global__ 
void surface_conservation_step(double * step, double const * vx, double const * vy, double const * vz, double const * vd, double const * lsf, double dt, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int rght1 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int rght2 	= sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int rght3 	= sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left1 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 	= sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3 	= sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

	double v1 = (lsf[left2] - lsf[left3]) / dx;
	double v2 = (lsf[left1] - lsf[left2]) / dx;
	double v3 = (lsf[ind]   - lsf[left1]) / dx;
	double v4 = (lsf[rght1] - lsf[ind])   / dx;
	double v5 = (lsf[rght2] - lsf[rght1]) / dx;
	double v6 = (lsf[rght3] - lsf[rght2]) / dx;

	double xL= weno_onesided_derivative(v1,v2,v3,v4,v5);
	double xR= weno_onesided_derivative(v6,v5,v4,v3,v2);

	int frnt1 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int frnt2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);	
	int frnt3 	= sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);	
	int back1 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3 	= sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

	v1 = (lsf[back2] - lsf[back3]) / dy;
	v2 = (lsf[back1] - lsf[back2]) / dy;
	v3 = (lsf[ind]	 - lsf[back1]) / dy;
	v4 = (lsf[frnt1] - lsf[ind])  / dy;
	v5 = (lsf[frnt2] - lsf[frnt1]) / dy;
	v6 = (lsf[frnt3] - lsf[frnt2]) / dy;
	
	double yB = weno_onesided_derivative(v1,v2,v3,v4,v5);
	double yF = weno_onesided_derivative(v6,v5,v4,v3,v2);

	int upup1	= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int upup2	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int upup3	= sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);	
	int down1 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 	= sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

	v1 = (lsf[down2] - lsf[down3]) / dz;
	v2 = (lsf[down1] - lsf[down2]) / dz;
	v3 = (lsf[ind]	 - lsf[down1]) / dz;
	v4 = (lsf[upup1] - lsf[ind])   / dz;
	v5 = (lsf[upup2] - lsf[upup1]) / dz;
	v6 = (lsf[upup3] - lsf[upup2]) / dz;

	double zD = weno_onesided_derivative(v1,v2,v3,v4,v5);
	double zU = weno_onesided_derivative(v6,v5,v4,v3,v2);

	step[ind] = (min2(0,vx[ind]) * xR + max2(0,vx[ind]) * xL + 
				 min2(0,vy[ind]) * yF + max2(0,vy[ind]) * yB + 
				 min2(0,vz[ind]) * zU + max2(0,vz[ind]) * zD +
				 lsf[ind] * vd[ind]) * dt ;
}	

// calculate numerical Hamiltonian for surface conservation law equation with finite volume method
// assuming c field and velocity field are already extended
__global__
void spatial_finite_volume_step(double * step, double const * vx, double const * vy, double const * vz, double const * lsf, double dt, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	double numericalFlux = 0; 
	double v_upwind, v_dnwind; // speed in the upwind and downwind direction
	int dnwind, upwind;

	// use linear approximation to calculate speed at the boundary
	dnwind 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	upwind 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);

	v_dnwind = (vx[dnwind] + vx[ind]) / 2.0;
	v_upwind = (vx[upwind] + vx[ind]) / 2.0; 

	numericalFlux += (min2(0,v_dnwind) * lsf[dnwind] - max2(0,v_dnwind) * lsf[ind]) * dy * dz;
	numericalFlux += (max2(0,v_upwind) * lsf[upwind] - min2(0,v_upwind) * lsf[ind]) * dy * dz;

	dnwind 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	upwind 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);

	v_dnwind = (vy[dnwind] + vy[ind]) / 2.0;
	v_upwind = (vy[upwind] + vy[ind]) / 2.0; 

	numericalFlux += (min2(0,v_dnwind) * lsf[dnwind] - max2(0,v_dnwind) * lsf[ind]) * dx * dz;
	numericalFlux += (max2(0,v_upwind) * lsf[upwind] - min2(0,v_upwind) * lsf[ind]) * dx * dz;

	dnwind	= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	upwind 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);

	v_dnwind = (vz[dnwind] + vz[ind]) / 2.0;
	v_upwind = (vz[upwind] + vz[ind]) / 2.0; 

	numericalFlux += (min2(0,v_dnwind) * lsf[dnwind] - max2(0,v_dnwind) * lsf[ind]) * dx * dy;
	numericalFlux += (max2(0,v_upwind) * lsf[upwind] - min2(0,v_upwind) * lsf[ind]) * dx * dy;

	step[ind] = numericalFlux * dt / (dx * dy * dz); 
}





















