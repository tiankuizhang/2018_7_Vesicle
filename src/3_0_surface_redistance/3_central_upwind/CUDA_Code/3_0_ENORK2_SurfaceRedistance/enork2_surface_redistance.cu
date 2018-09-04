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
void advection_velocity(double & H1, double & H2, double & H3, double & normal_d, double sign, double Dx, double Dy, double Dz, double nx, double ny, double nz)
{
	normal_d = nx * Dx + ny * Dy + nz * Dz;

	H1 = sign * (Dx - nx * normal_d);
	H2 = sign * (Dy - ny * normal_d);
	H3 = sign * (Dz - nz * normal_d);

	double H_mag = sqrt(H1*H1+H2*H2+H3*H3)+1e-10;

	H1 = H1/H_mag;
	H2 = H2/H_mag;
	H3 = H3/H_mag;
}

__device__ inline
void Upwind_Hamiltonian(double & Hamil, double normal_d, double sign, double Dx, double Dy, double Dz)
{
	Hamil = sign* ( sqrt( Dx*Dx + Dy*Dy + Dz*Dz - normal_d*normal_d ) - 1);
	//Hamil =   Dx*Dx + Dy*Dy + Dz*Dz - normal_d*normal_d  ;
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

// calculate surface redistance step with central upwind scheme
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

	double_eno_derivative eno_dy = eno_derivative( lsf[back2], lsf[back], lsf[ind], lsf[front], lsf[front2], ypf[ind], ypb[ind], dy);
	double Dy[2] = {eno_dy.sR, eno_dy.sL};

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
	//double epsilon=1e-10;
	double epsilon=0.;
	double numerical_Hamiltonian = 0;
	double denominator = (a[0]+a[1])*(b[0]+b[1])*(c[0]+c[1])+epsilon;

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

	numerical_Hamiltonian += - a[0]*a[1]*(Dx[0]-Dx[1])/(a[0]+a[1]+epsilon) - b[0]*b[1]*(Dy[0]-Dy[1])/(b[0]+b[1]+epsilon) - c[0]*c[1]*(Dz[0]-Dz[1])/(c[0]+c[1]+epsilon);

	step[ind] = numerical_Hamiltonian * deltat[ind];

	//step[ind] = numerical_Hamiltonian;
	//step[ind] = Hamiltonian[0];

}









