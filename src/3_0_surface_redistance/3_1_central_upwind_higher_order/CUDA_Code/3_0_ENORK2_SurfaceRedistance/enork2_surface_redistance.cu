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
void advection_velocity(double & H1, double & H2, double & H3, double & normal_d, double sign, double Dx, double Dy, double Dz, double nx, double ny, double nz)
{
	normal_d = nx * Dx + ny * Dy + nz * Dz;

	H1 = sign * (Dx - nx * normal_d);
	H2 = sign * (Dy - ny * normal_d);
	H3 = sign * (Dz - nz * normal_d);

	double H_mag = sqrt(H1*H1+H2*H2+H3*H3);

	H1 = H1/H_mag;
	H2 = H2/H_mag;
	H3 = H3/H_mag;
}

__device__ inline
void Upwind_Hamiltonian(double & Hamil, double normal_d, double sign, double Dx, double Dy, double Dz)
{
	// numerical error can lead to negative value inside sqrt() 
	// the following code is needed to avoid NAN due to sqrt of a negative number
	Hamil = sign* ( sqrt( max2(0,Dx*Dx+Dy*Dy+Dz*Dz-normal_d*normal_d) ) - 1);
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

/****************************************************************************** 
 * calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
 * without boundary correction
 ******************************************************************************/
__device__ inline
double_eno_derivative eno_derivative_field( double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m;
	double_eno_derivative eno_d;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	eno_d.sR = (v2 - v0) / ds - ds * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	eno_d.sL = (v0 - v1) / ds + ds * p2m;

	return eno_d;

}

// calculate surface redistance step with central upwind scheme without higher order term
// now lsf represents the auxilary level set function(not the level set function)
// inputs : the auxilary level set function, sign of the initial level set function, distance to the interface, normal vectors
__global__
void surface_redistance_step_1st(double * step, double const * lsf, double const * sign, double const * deltat, double const * nx, double const * ny, double const * nz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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
	double Dx[2] = {eno_dx.sR, eno_dx.sL};

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy); double Dy[2] = {eno_dy.sR, eno_dy.sL};

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	double Dz[2] = {eno_dz.sR, eno_dz.sL};

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];
	double Sign = sign[ind];

	// different choices yield different upwind direction
	int const choice_x[8]={0,0,0,0,1,1,1,1};
	int const choice_y[8]={0,0,1,1,0,0,1,1};
	int const choice_z[8]={0,1,0,1,0,1,0,1};

	double Hamiltonian[8]={0,0,0,0,0,0,0,0};

	// a[0],a[1] is the magnitude of the maximum forward and backward
    // information propagation speed in the x direction. b for y and c for z direction
	double a[2]={0,0};
	double b[2]={0,0};
	double c[2]={0,0};

	for(int i=0;i<8;i++){
		double dr_x = Dx[choice_x[i]];
		double dr_y = Dy[choice_y[i]];
		double dr_z = Dz[choice_z[i]];

		double H1,H2,H3, normal_d;
		advection_velocity(H1,H2,H3,normal_d,Sign,dr_x,dr_y,dr_z,Nx,Ny,Nz);
		Upwind_Hamiltonian(Hamiltonian[i],normal_d,Sign,dr_x,dr_y,dr_z);	
		
		a[0] = max2(a[0],max2(H1,0));
		b[0] = max2(b[0],max2(H2,0));
		c[0] = max2(c[0],max2(H3,0));

		a[1] = fabs(min2(-a[1],min2(H1,0)));
		b[1] = fabs(min2(-b[1],min2(H2,0)));
		c[1] = fabs(min2(-c[1],min2(H3,0)));
	}

	// calculate the numerical Hamiltonian
	//double epsilon=1e-6;
	double numerical_Hamiltonian = 0;
	double denominator = (a[0]+a[1])*(b[0]+b[1])*(c[0]+c[1]);

	int const choice_a[8]={1,1,1,1,0,0,0,0};
	int const choice_b[8]={1,1,0,0,1,1,0,0};
	int const choice_c[8]={1,0,1,0,1,0,1,0};

	for(int i=0;i<8;i++){
		double	H_a = a[choice_a[i]];
		double	H_b = b[choice_b[i]];
		double	H_c = c[choice_c[i]];
		numerical_Hamiltonian += H_a * H_b * H_c * Hamiltonian[i];
	}
	numerical_Hamiltonian = numerical_Hamiltonian/denominator;

	numerical_Hamiltonian += - a[0]*a[1]*(Dx[0]-Dx[1])/(a[0]+a[1]) - b[0]*b[1]*(Dy[0]-Dy[1])/(b[0]+b[1]) - c[0]*c[1]*(Dz[0]-Dz[1])/(c[0]+c[1]);

	step[ind] = numerical_Hamiltonian * deltat[ind];

}



// calculate surface redistance step with central upwind scheme with higher order term
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
	double Dx[2] = {eno_dx.sR, eno_dx.sL};

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy); double Dy[2] = {eno_dy.sR, eno_dy.sL};

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	double Dz[2] = {eno_dz.sR, eno_dz.sL};

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];
	double Sign = sign[ind];

	// different choices yield different upwind direction
	int const choice_x[8]={0,0,0,0,1,1,1,1};
	int const choice_y[8]={0,0,1,1,0,0,1,1};
	int const choice_z[8]={0,1,0,1,0,1,0,1};

	double Hamiltonian[8]={0,0,0,0,0,0,0,0};

	// a[0],a[1] is the magnitude of the maximum forward and backward
    // information propagation speed in the x direction. b for y and c for z direction
	double a[2]={0,0};
	double b[2]={0,0};
	double c[2]={0,0};

	for(int i=0;i<8;i++){
		double dr_x = Dx[choice_x[i]];
		double dr_y = Dy[choice_y[i]];
		double dr_z = Dz[choice_z[i]];

		double H1,H2,H3, normal_d;
		advection_velocity(H1,H2,H3,normal_d,Sign,dr_x,dr_y,dr_z,Nx,Ny,Nz);
		Upwind_Hamiltonian(Hamiltonian[i],normal_d,Sign,dr_x,dr_y,dr_z);	
		
		a[0] = max2(a[0],max2(H1,0));
		b[0] = max2(b[0],max2(H2,0));
		c[0] = max2(c[0],max2(H3,0));

		a[1] = fabs(min2(-a[1],min2(H1,0)));
		b[1] = fabs(min2(-b[1],min2(H2,0)));
		c[1] = fabs(min2(-c[1],min2(H3,0)));
	}

	// calculate the numerical Hamiltonian
	//double epsilon=1e-6;
	double numerical_Hamiltonian = 0;
	double denominator = (a[0]+a[1])*(b[0]+b[1])*(c[0]+c[1]);

	int const choice_a[8]={1,1,1,1,0,0,0,0};
	int const choice_b[8]={1,1,0,0,1,1,0,0};
	int const choice_c[8]={1,0,1,0,1,0,1,0};

	for(int i=0;i<8;i++){
		double	H_a = a[choice_a[i]];
		double	H_b = b[choice_b[i]];
		double	H_c = c[choice_c[i]];
		numerical_Hamiltonian += H_a * H_b * H_c * Hamiltonian[i];
	}
	numerical_Hamiltonian = numerical_Hamiltonian/denominator;

	numerical_Hamiltonian += - a[0]*a[1]*(Dx[0]-Dx[1])/(a[0]+a[1]) - b[0]*b[1]*(Dy[0]-Dy[1])/(b[0]+b[1]) - c[0]*c[1]*(Dz[0]-Dz[1])/(c[0]+c[1]);

	//
	// calculate higher order terms
	double psi_x_int_pp = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[0]-Hamiltonian[4]) ) / (a[0]+a[1]) ;
	double psi_x_int_mm = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[3]-Hamiltonian[7]) ) / (a[0]+a[1]) ;
	double psi_x_int_pm = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[1]-Hamiltonian[5]) ) / (a[0]+a[1]) ;
	double psi_x_int_mp = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[2]-Hamiltonian[6]) ) / (a[0]+a[1]) ;

	double dx_pp = min_mod(Dx[0]-psi_x_int_pp,psi_x_int_pp-Dx[1]);
	double dx_mm = min_mod(Dx[0]-psi_x_int_mm,psi_x_int_mm-Dx[1]);
	double dx_pm = min_mod(Dx[0]-psi_x_int_pm,psi_x_int_pm-Dx[1]);
	double dx_mp = min_mod(Dx[0]-psi_x_int_mp,psi_x_int_mp-Dx[1]);

	double x_c = b[0]*c[0]*dx_mm + b[1]*c[1]*dx_pp + b[0]*c[1]*dx_mp + b[1]*c[0]*dx_pm;

	double psi_y_int_pp = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[0]-Hamiltonian[2]) ) / (b[0]+b[1]) ;
	double psi_y_int_mm = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[5]-Hamiltonian[7]) ) / (b[0]+b[1]) ;
	double psi_y_int_pm = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[1]-Hamiltonian[3]) ) / (b[0]+b[1]) ;
	double psi_y_int_mp = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[4]-Hamiltonian[6]) ) / (b[0]+b[1]) ;

	double dy_pp = min_mod(Dy[0]-psi_y_int_pp,psi_y_int_pp-Dy[1]);
	double dy_mm = min_mod(Dy[0]-psi_y_int_mm,psi_y_int_mm-Dy[1]);
	double dy_pm = min_mod(Dy[0]-psi_y_int_pm,psi_y_int_pm-Dy[1]);
	double dy_mp = min_mod(Dy[0]-psi_y_int_mp,psi_y_int_mp-Dy[1]);

	double y_c = a[0]*c[0]*dy_mm + a[1]*c[1]*dy_pp + a[0]*c[1]*dy_mp + a[1]*c[0]*dy_pm;

	double psi_z_int_pp = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[0]-Hamiltonian[1]) ) / (c[0]+c[1]) ;
	double psi_z_int_mm = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[6]-Hamiltonian[7]) ) / (c[0]+c[1]) ;
	double psi_z_int_pm = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[2]-Hamiltonian[3]) ) / (c[0]+c[1]) ;
	double psi_z_int_mp = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[4]-Hamiltonian[5]) ) / (c[0]+c[1]) ;

	double dz_pp = min_mod(Dz[0]-psi_z_int_pp,psi_z_int_pp-Dz[1]);
	double dz_mm = min_mod(Dz[0]-psi_z_int_mm,psi_z_int_mm-Dz[1]);
	double dz_pm = min_mod(Dz[0]-psi_z_int_pm,psi_z_int_pm-Dz[1]);
	double dz_mp = min_mod(Dz[0]-psi_z_int_mp,psi_z_int_mp-Dz[1]);

	double z_c = a[0]*b[0]*dz_mm + a[1]*b[1]*dz_pp + a[0]*b[1]*dz_mp + a[1]*b[0]*dz_pm;

	numerical_Hamiltonian += (a[0]*a[1]*x_c + b[0]*b[1]*y_c + c[0]*c[1]*z_c) / denominator;
	//

	step[ind] = numerical_Hamiltonian * deltat[ind];

}


// calculate surface redistance step with central upwind scheme with higher order term
// AND WITHOUT boundary correction
// now lsf represents the auxilary level set function(not the level set function)
// inputs : the auxilary level set function, sign of the initial level set function, distance to the interface, normal vectors
__global__
void surface_redistance_step_noboundaryfix(double * step, double const * lsf, double const * sign, double const * deltat, double const * nx, double const * ny, double const * nz, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
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

	double_eno_derivative eno_dx = eno_derivative_field( lsf[left2], lsf[left], lsf[ind], lsf[right], lsf[right2], xpr[ind], xpl[ind], dx);
	double Dx[2] = {eno_dx.sR, eno_dx.sL};

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	double_eno_derivative eno_dy = eno_derivative_field( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy); double Dy[2] = {eno_dy.sR, eno_dy.sL};

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double_eno_derivative eno_dz = eno_derivative_field( lsf[down2], lsf[down], lsf[ind], lsf[up], lsf[up2], zpu[ind], zpd[ind], dz);
	double Dz[2] = {eno_dz.sR, eno_dz.sL};

	double Nx = nx[ind];
	double Ny = ny[ind];
	double Nz = nz[ind];
	double Sign = sign[ind];

	// different choices yield different upwind direction
	int const choice_x[8]={0,0,0,0,1,1,1,1};
	int const choice_y[8]={0,0,1,1,0,0,1,1};
	int const choice_z[8]={0,1,0,1,0,1,0,1};

	double Hamiltonian[8]={0,0,0,0,0,0,0,0};

	// a[0],a[1] is the magnitude of the maximum forward and backward
    // information propagation speed in the x direction. b for y and c for z direction
	double a[2]={0,0};
	double b[2]={0,0};
	double c[2]={0,0};

	for(int i=0;i<8;i++){
		double dr_x = Dx[choice_x[i]];
		double dr_y = Dy[choice_y[i]];
		double dr_z = Dz[choice_z[i]];

		double H1,H2,H3, normal_d;
		advection_velocity(H1,H2,H3,normal_d,Sign,dr_x,dr_y,dr_z,Nx,Ny,Nz);
		Upwind_Hamiltonian(Hamiltonian[i],normal_d,Sign,dr_x,dr_y,dr_z);	
		
		a[0] = max2(a[0],max2(H1,0));
		b[0] = max2(b[0],max2(H2,0));
		c[0] = max2(c[0],max2(H3,0));

		a[1] = fabs(min2(-a[1],min2(H1,0)));
		b[1] = fabs(min2(-b[1],min2(H2,0)));
		c[1] = fabs(min2(-c[1],min2(H3,0)));
	}

	// calculate the numerical Hamiltonian
	//double epsilon=1e-6;
	double numerical_Hamiltonian = 0;
	double denominator = (a[0]+a[1])*(b[0]+b[1])*(c[0]+c[1]);

	int const choice_a[8]={1,1,1,1,0,0,0,0};
	int const choice_b[8]={1,1,0,0,1,1,0,0};
	int const choice_c[8]={1,0,1,0,1,0,1,0};

	for(int i=0;i<8;i++){
		double	H_a = a[choice_a[i]];
		double	H_b = b[choice_b[i]];
		double	H_c = c[choice_c[i]];
		numerical_Hamiltonian += H_a * H_b * H_c * Hamiltonian[i];
	}
	numerical_Hamiltonian = numerical_Hamiltonian/denominator;

	numerical_Hamiltonian += - a[0]*a[1]*(Dx[0]-Dx[1])/(a[0]+a[1]) - b[0]*b[1]*(Dy[0]-Dy[1])/(b[0]+b[1]) - c[0]*c[1]*(Dz[0]-Dz[1])/(c[0]+c[1]);

	//
	// calculate higher order terms
	double psi_x_int_pp = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[0]-Hamiltonian[4]) ) / (a[0]+a[1]) ;
	double psi_x_int_mm = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[3]-Hamiltonian[7]) ) / (a[0]+a[1]) ;
	double psi_x_int_pm = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[1]-Hamiltonian[5]) ) / (a[0]+a[1]) ;
	double psi_x_int_mp = ( (a[0]*Dx[0]+a[1]*Dx[1]) - (Hamiltonian[2]-Hamiltonian[6]) ) / (a[0]+a[1]) ;

	double dx_pp = min_mod(Dx[0]-psi_x_int_pp,psi_x_int_pp-Dx[1]);
	double dx_mm = min_mod(Dx[0]-psi_x_int_mm,psi_x_int_mm-Dx[1]);
	double dx_pm = min_mod(Dx[0]-psi_x_int_pm,psi_x_int_pm-Dx[1]);
	double dx_mp = min_mod(Dx[0]-psi_x_int_mp,psi_x_int_mp-Dx[1]);

	double x_c = b[0]*c[0]*dx_mm + b[1]*c[1]*dx_pp + b[0]*c[1]*dx_mp + b[1]*c[0]*dx_pm;

	double psi_y_int_pp = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[0]-Hamiltonian[2]) ) / (b[0]+b[1]) ;
	double psi_y_int_mm = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[5]-Hamiltonian[7]) ) / (b[0]+b[1]) ;
	double psi_y_int_pm = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[1]-Hamiltonian[3]) ) / (b[0]+b[1]) ;
	double psi_y_int_mp = ( (b[0]*Dy[0]+b[1]*Dy[1]) - (Hamiltonian[4]-Hamiltonian[6]) ) / (b[0]+b[1]) ;

	double dy_pp = min_mod(Dy[0]-psi_y_int_pp,psi_y_int_pp-Dy[1]);
	double dy_mm = min_mod(Dy[0]-psi_y_int_mm,psi_y_int_mm-Dy[1]);
	double dy_pm = min_mod(Dy[0]-psi_y_int_pm,psi_y_int_pm-Dy[1]);
	double dy_mp = min_mod(Dy[0]-psi_y_int_mp,psi_y_int_mp-Dy[1]);

	double y_c = a[0]*c[0]*dy_mm + a[1]*c[1]*dy_pp + a[0]*c[1]*dy_mp + a[1]*c[0]*dy_pm;

	double psi_z_int_pp = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[0]-Hamiltonian[1]) ) / (c[0]+c[1]) ;
	double psi_z_int_mm = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[6]-Hamiltonian[7]) ) / (c[0]+c[1]) ;
	double psi_z_int_pm = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[2]-Hamiltonian[3]) ) / (c[0]+c[1]) ;
	double psi_z_int_mp = ( (c[0]*Dz[0]+c[1]*Dz[1]) - (Hamiltonian[4]-Hamiltonian[5]) ) / (c[0]+c[1]) ;

	double dz_pp = min_mod(Dz[0]-psi_z_int_pp,psi_z_int_pp-Dz[1]);
	double dz_mm = min_mod(Dz[0]-psi_z_int_mm,psi_z_int_mm-Dz[1]);
	double dz_pm = min_mod(Dz[0]-psi_z_int_pm,psi_z_int_pm-Dz[1]);
	double dz_mp = min_mod(Dz[0]-psi_z_int_mp,psi_z_int_mp-Dz[1]);

	double z_c = a[0]*b[0]*dz_mm + a[1]*b[1]*dz_pp + a[0]*b[1]*dz_mp + a[1]*b[0]*dz_pm;

	numerical_Hamiltonian += (a[0]*a[1]*x_c + b[0]*b[1]*y_c + c[0]*c[1]*z_c) / denominator;
	//

	step[ind] = numerical_Hamiltonian * deltat[ind];

}















