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
double max_abs(double x, double y)
{
	return (fabs(x)<fabs(y) ? y : x);
}

__device__ inline
double sign(double x)
{
	return (x>0) ? 1.0 : -1.0;
}

__device__ inline
void upwind_derivative(double & Dr, double Cx, double Cy, double Cz, double xR, double xL, double yF, double yB, double zU, double zD)
{
	Dr = 0;
	Dr += max2(Cx,0)*xR - min2(Cx,0)*xL;
	Dr += max2(Cy,0)*yF - min2(Cy,0)*yB;
	Dr += max2(Cz,0)*zU - min2(Cz,0)*zD;
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
void surface_redistance_step(double * step, double const * lsf, double const * sign, double const * deltat, double const * nx, double const * ny, double const * nz, double const * ax, double const * ay, double const * az, double const * bx, double const * by, double const * bz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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
	//double Dx[2] = {eno_dx.sR, eno_dx.sL};
	double xR = eno_dx.sR;
	double xL = eno_dx.sL;

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);
	//double Dy[2] = {eno_dy.sR, eno_dy.sL};
	double yF = eno_dy.sR;
	double yB = eno_dy.sL;

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	//double Dz[2] = {eno_dz.sR, eno_dz.sL};
	double zU = eno_dz.sR;
	double zD = eno_dz.sL;

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];
	double Ax = ax[ind];
	double Ay = ay[ind];
	double Az = az[ind];
	double Bx = bx[ind];
	double By = by[ind];
	double Bz = bz[ind];
	
	double Sign = sign[ind];
	double dt = deltat[ind];

	// forward and backward derivatives along a/b directions
	double Da_f, Da_b, Db_f, Db_b;
	upwind_derivative(Da_f,Ax,Ay,Az,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Da_b,-Ax,-Ay,-Az,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Db_f,Bx,By,Bz,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Db_b,-Bx,-By,-Bz,xR,xL,yF,yB,zU,zD);

	// choice of one-sided derivatives in the surface coordiante
	double Da = max_abs(max2(Sign*Da_b,0),min2(Sign*Da_f,0));
	double Db = max_abs(max2(Sign*Db_b,0),min2(Sign*Db_f,0));
	double den = sqrt(Da*Da+Db*Db);

	// direction of information 
	double H1 = (Da * Ax + Db * Bx) / den;
	double H2 = (Da * Ay + Db * By) / den;
	double H3 = (Da * Az + Db * Bz) / den;

	// select the right derivatives 
	double dr_x = (H1>0) ? xL : ( (H1<0) ? xR : 0 );
	double dr_y = (H2>0) ? yB : ( (H2<0) ? yF : 0 );
	double dr_z = (H3>0) ? zD : ( (H3<0) ? zU : 0 );

	double dr_n = Nx * dr_x + Ny * dr_y + Nz * dr_z;	
	step[ind] = dt * Sign * (sqrt( dr_x*dr_x + dr_y*dr_y + dr_z*dr_z - dr_n*dr_n) - 1);

}

// construct a local coordinate system at each point
__global__
void surface_coordinate(double * ax, double * ay, double * az, double * bx, double * by, double * bz, double const * nx, double const * ny, double const * nz, int rows, int cols, int pges)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];

	double epsilon =  1e-6;

	if(fabs(Nx)>epsilon){
		double den1 = sqrt(Nx*Nx+Ny*Ny);
		ax[ind] = -Ny/den1;
		ay[ind] = Nx/den1;
		
		den1 = sqrt(Nx*Nx+Nz*Nz);
		bx[ind] = - Nz/den1;
		bz[ind] = Nx/den1;
	}
	else if(fabs(Ny)>epsilon){
		double den1 = sqrt(Ny*Ny+Nx*Nx);
		ay[ind] = -Nx/den1;
		ax[ind] = Ny/den1;

		den1 = sqrt(Ny*Ny+Nz*Nz);
		by[ind] = - Nz/den1;
	    bz[ind] = Ny/den1;	
	}
	else if(fabs(Nz)>epsilon){
		double den1 = sqrt(Nz*Nz+Nx*Nx);
		az[ind] = -Nx/den1;
		ax[ind] = Nz/den1;
		
		den1 = sqrt(Nz*Nz+Ny*Ny);
		bz[ind] = - Ny/den1;
		by[ind] = Nz/den1;
	}


}

__global__
void surface_redistance_step_backup(double * step, double const * lsf, double const * sign, double const * deltat, double const * ax, double const * ay, double const * az, double const * bx, double const * by, double const * bz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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
	//double Dx[2] = {eno_dx.sR, eno_dx.sL};
	double xR = eno_dx.sR;
	double xL = eno_dx.sL;

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);
	//double Dy[2] = {eno_dy.sR, eno_dy.sL};
	double yF = eno_dy.sR;
	double yB = eno_dy.sL;

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	//double Dz[2] = {eno_dz.sR, eno_dz.sL};
	double zU = eno_dz.sR;
	double zD = eno_dz.sL;

	double Ax = ax[ind];
	double Ay = ay[ind];
	double Az = az[ind];
	double Bx = bx[ind];
	double By = by[ind];
	double Bz = bz[ind];
	
	double Sign = sign[ind];
	double dt = deltat[ind];

	// forward and backward derivatives along a/b directions
	double Da_f, Da_b, Db_f, Db_b;
	upwind_derivative(Da_f,Ax,Ay,Az,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Da_b,-Ax,-Ay,-Az,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Db_f,Bx,By,Bz,xR,xL,yF,yB,zU,zD);
	upwind_derivative(Db_b,-Bx,-By,-Bz,xR,xL,yF,yB,zU,zD);

	double Da = max(max2(Sign*Da_b,0),-min2(Sign*Da_f,0));
	double Db = max(max2(Sign*Db_b,0),-min2(Sign*Db_f,0));

	step[ind] = dt * Sign * (sqrt(Da*Da + Db*Db) - 1);

}








