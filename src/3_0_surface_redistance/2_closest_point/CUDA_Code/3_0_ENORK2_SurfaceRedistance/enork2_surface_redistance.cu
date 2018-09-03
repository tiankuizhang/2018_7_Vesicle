/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the surface redistance scheme 
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

__device__ inline
bool same_sign(double x, double y)
{
	return (x*y>0) || (x==0 && y==0);
}	

__device__ inline
void advection_velocity(double & H1, double & H2, double & H3, double sign, double Dx, double Dy, double Dz, double nx, double ny, double nz)
{
	double normal_d = nx * Dx + ny + Dy + nz * Dz;

	H1 = sign * (Dx - nx * normal_d);
	H2 = sign * (Dy - ny * normal_d);
	H3 = sign * (Dz - nz * normal_d);

	double H_mag = sqrt(H1*H1+H2*H2+H3*H3+1e-6);

	H1 = H1/H_mag;
	H2 = H2/H_mag;
	H3 = H3/H_mag;

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
 ******************************************************************************/
__device__ inline
double_eno_derivative eno_derivative( double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m;
	double_eno_derivative eno_d;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	double vr = (pr==ds) ? v2 : 0;
	eno_d.sR = (vr - v0) / pr - pr * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	double vl = (pl==ds) ? v1 : 0;
	eno_d.sL = (v0 - vl) / pl + pl * p2m;

	return eno_d;

}

// calculate surface redistance step
// now lsf represents the auxilary level set function(not the level set function)
// inputs : the auxilary level set function, sign of the initial level set function, distance to the interface, normal vectors
__global__
void surface_redistance_step(double * step, double const * lsf, double const * sign, double const * deltat, double const * nx, double const * ny, double const * nz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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

	double_eno_derivative eno_dx = eno_derivative( lsf[left2], lsf[left], lsf[ind], lsf[right], lsf[right2], xpr[ind], xpl[ind], dx);
	double Dx[3] = {eno_dx.sR, 0, eno_dx.sL};

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);
	double Dy[3] = {eno_dy.sR, 0, eno_dy.sL};

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	double Dz[3] = {eno_dz.sR, 0, eno_dz.sL};


	//Forward=-1, None=0, BackWard=1
	int const choice_x[26] = {-1,-1,-1,-1, 1, 1, 1, 1,	 0, 0, 0, 0,-1,-1, 1, 1,-1,-1, 1, 1,	-1, 1, 0, 0, 0, 0};
	int const choice_y[26] = {-1,-1, 1, 1,-1,-1, 1, 1,	-1,-1, 1, 1, 0, 0, 0, 0,-1, 1,-1, 1,	 0, 0,-1, 1, 0, 0};
	int const choice_z[26] = {-1, 1,-1, 1,-1, 1,-1, 1,	-1, 1,-1, 1,-1, 1,-1, 1, 0, 0, 0, 0,	 0, 0, 0, 0,-1, 1};

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];
	double Sign = sign[ind];

	double dx_c = (Dx[0] + Dx[2]) / 2;
	double dy_c = (Dy[0] + Dy[2]) / 2;
	double dz_c = (Dz[0] + Dz[2]) / 2;

	double maxH1 = 0;
	double maxH2 = 0;
	double maxH3 = 0;

	for(int i=0;i<26;i++){
		double dr_x = Dx[choice_x[i]+1];
		double dr_y = Dy[choice_y[i]+1];
		double dr_z = Dz[choice_z[i]+1];

		double H1, H2, H3; // information propagation direction
		advection_velocity(H1,H2,H3,Sign,dr_x,dr_y,dr_z,Nx,Ny,Nz);

		maxH1 = (fabs(H1)>maxH1) ? fabs(H1) : maxH1;
		maxH2 = (fabs(H2)>maxH2) ? fabs(H2) : maxH1;
		maxH3 = (fabs(H3)>maxH3) ? fabs(H3) : maxH1;

	}


	double dt = deltat[ind];
	//step[ind] = dt*Sign*(sqrt( pow(dx_c*Nz-Nx*dz_c,2)+pow(dy_c*Nx-Ny*dx_c,2)+pow(dz_c*Ny-Nz*dy_c,2) )-1) - 0.5*(Dx[0]-Dx[2]) - 0.5*(Dy[0]-Dy[2]) - 0.5*(Dz[0]-Dz[2]);
	step[ind] = dt*Sign*(sqrt( pow(dx_c*Nz-Nx*dz_c,2)+pow(dy_c*Nx-Ny*dx_c,2)+pow(dz_c*Ny-Nz*dy_c,2) )-1) - 0.5*dt*(maxH1*(Dx[0]-Dx[2]) - maxH2*(Dy[0]-Dy[2]) - maxH3*(Dz[0]-Dz[2]));
	//step[ind] = deltat[ind]*Sign*(sqrt( pow(dx_c*Nz-Nx*dz_c,2)+pow(dy_c*Nx-Ny*dx_c,2)+pow(dz_c*Ny-Nz*dy_c,2) )-1);
	//step[ind] = maxH3;
	//step[ind] = Dz[0] - Dz[2];
//	step[ind] = Dz[0];


}









