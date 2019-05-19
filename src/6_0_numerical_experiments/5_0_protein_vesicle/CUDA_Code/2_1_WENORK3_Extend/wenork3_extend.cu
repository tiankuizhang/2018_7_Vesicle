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
double min_abs(double x, double y)
{
	return (fabs(x)<fabs(y)) ? x : y;
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

// given a stencil across the boundary: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// create a new stencil (x3m,h3m),(x2m,h2m),(x1m,h1m),(x0,h0),(x1,h1),(x2,h2),(x3,h3) including boundary nodes
__device__ inline
void select_stencil(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds, double h_fore, double h_back)
{
		h0 = p4; x0 = 0.0;

		if(r1<ds){
			x1 = r1;
			x2 = ds;
			x3 = 2*ds;
			h1 = h_fore;
			h2 = p5;
			h3 = p6;
		}else{
			x1 = ds;
			x2 = 2*ds;
			x3 = 3*ds;
			h1 = p5;
			h2 = p6;
			h3 = p7;
		}

		if(l1<ds){
			x1m = -l1;
			x2m = - ds;
			x3m = - 2*ds;
			h1m = h_back;
			h2m = p3;
			h3m = p2;
		}else{
			x1m = -ds;
			x2m = - 2*ds;
			x3m = - 3*ds;
			h1m = p3;
			h2m = p2;
			h3m = p1;
		}
}

// for stencil (x3m,h3m),(x2m,h2m),(x1m,h1m),(x0,h0),(x1,h1),(x2,h2),(x3,h3) including boundary nodes
// calculate cubic eno derivatives at (x0,h0)
// note that this is a nonuniform stencil
__device__ inline
void ENO_cubic_derivative(double & d_fore, double & d_back, double h3m, double h2m, double h1m, double h0, double h1, double h2, double h3, double x3m, double x2m, double x1m, double x0, double x1, double x2, double x3)
{

	// divided differences
	double d1_2_5   = (h3  - h2)  / (x3  - x2) ;
	double d1_1_5   = (h2  - h1)  / (x2  - x1) ;
	double d1_0_5   = (h1  - h0)  / (x1  - x0) ;
	double d1_m0_5  = (h0  - h1m) / (x0  - x1m);
	double d1_m1_5  = (h1m - h2m) / (x1m - x2m);
	double d1_m2_5  = (h2m - h3m) / (x2m - x3m);

	double d2_2  = (d1_2_5  - d1_1_5)  / (x3  - x1) ;
	double d2_1  = (d1_1_5  - d1_0_5)  / (x2  - x0) ;
	double d2_0  = (d1_0_5  - d1_m0_5) / (x1  - x1m);
	double d2_m1 = (d1_m0_5 - d1_m1_5) / (x0  - x2m);
	double d2_m2 = (d1_m1_5 - d1_m2_5) / (x1m - x3m);

	double d3_1_5  = (d2_2  - d2_1)  / (x3 - x0) ;
	double d3_0_5  = (d2_1  - d2_0)  / (x2 - x1m);
	double d3_m0_5 = (d2_0  - d2_m1) / (x1 - x2m);
	double d3_m1_5 = (d2_m1 - d2_m2) / (x0 - x3m);

	double a1 = (x0 - x1m) * (x0 - x2m) * min_abs(d3_m0_5, d3_m1_5);
	double a2 = (x0 - x1m) * (x0 - x1)  * min_abs(d3_m0_5, d3_0_5);
	double a = (fabs(d2_m1) < fabs(d2_0)) ? a1 : a2;

	double b1 = (x0 - x1m) * (x0 - x1) * min_abs(d3_m0_5, d3_0_5);
	double b2 = (x0 - x1)  * (x0 - x2) * min_abs(d3_0_5,  d3_1_5);
	double b = (fabs(d2_0) < fabs(d2_1)) ? b1 : b2;

	d_back = d1_m0_5 + min_mod(d2_m1,d2_0) * (x0 - x1m) + a;
	d_fore = d1_0_5  + min_mod(d2_0, d2_1) * (x0 - x1)  + b;
}

// calculate weno derivative at p4: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// where px are level set function values at node x
// lx, rx are distance to the left/right node 
__device__ inline
void weno_derivative_boundary(double & d_fore, double & d_back, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds, double v_fore, double v_back)
{
	// the condistion below is better than p3*p4<0 || p4*p5<0 
	bool cross_interface = r1<ds || l1<ds;

	if(!cross_interface){
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p5) / ds;
		double v6 = (p7 - p6) / ds;
		d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
	}// if not a node IMMEDIATELY adjacent to the boundary, calculate weno derivatives as usual
	else{
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds,v_fore,v_back);
		ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
	}// for nodes IMMEDIATELY adjacent to the boundary, use cubic ENO interpolant
}

/********************************************************************************
********************************************************************************/
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
/********************************************************************************
********************************************************************************/

__device__ inline
void cubic_interp_coefficient(double & c0, double & c1, double & c2, double & c3, double v0, double v1, double v2, double v3, double s)
{
	c0 = v0;
 	c1 = (     3.0 * (v1-v0) - 3.0/2.0 * (v2-v0) + 1.0/3.0 * (v3-v0) ) / s;
 	c2 = (-5.0/2.0 * (v1-v0) +     2.0 * (v2-v0) - 1.0/2.0 * (v3-v0) ) / pow(s,2);
 	c3 = ( 1.0/2.0 * (v1-v0) - 1.0/2.0 * (v2-v0) + 1.0/6.0 * (v3-v0) ) / pow(s,3);
}

// modify forward/backward(c_f/b) values at v0:[v4,v1,v0,v2,v3]
// if dis_b/f!=ds, then there is boundary nearby, a cubic interpolant is then constructed
// C(x) = c0 + c1*x + c2*x^2 + c3*c^3 through (0,v1),(ds,v0),(2*ds,v2),(3*ds,v2) assuming boundary is between v0,v2
// and used to calculate c_b/f at boundary crossing nodes
__device__ inline
void cubic_interp(double & c_forward, double & c_backward, double dis_f, double dis_b, double ds, double v4, double v1, double v0, double v2, double v3)
{
	c_forward = 0.;
	c_backward = 0.;

	double c0,c1,c2,c3; // coefficient for cubic interpolant
	// if there is a boundary in the forward direction
	if(dis_f<ds){
		cubic_interp_coefficient(c0,c1,c2,c3,v1,v0,v2,v3,ds);
		double xc = ds + dis_f; // coordinate of the boundary point
		c_forward = c0 + c1 * xc + c2 * pow(xc,2) + c3 * pow(xc,3);
	}
	// if there is a boundary in the backward direction
	if(dis_b<ds){
		cubic_interp_coefficient(c0,c1,c2,c3,v4,v1,v0,v2,ds);
		double xc = 2*ds - dis_b; 
		c_backward = c0 + c1 * xc + c2 * pow(xc,2) + c3 * pow(xc,3);
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

	cubic_interp(cpr[ind],cpl[ind],xpr[ind],xpl[ind],dx,lsf[left2],lsf[left],lsf[ind],lsf[right],lsf[right2]);

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int front2 	= sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 	= sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	cubic_interp(cpf[ind],cpb[ind],ypf[ind],ypb[ind],dy,lsf[back2],lsf[back],lsf[ind],lsf[front],lsf[front2]);

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int up2 	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	cubic_interp(cpu[ind],cpd[ind],zpu[ind],zpd[ind],dz,lsf[down2],lsf[down],lsf[ind],lsf[up],lsf[up2]);

}
/********************************************************************************
********************************************************************************/


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

	double epsilon = 1e-6 * dx;
	if( xpr[ind]<epsilon || xpl[ind]<epsilon ||
		ypf[ind]<epsilon || ypb[ind]<epsilon || 
	   	zpu[ind]<epsilon || zpd[ind]<epsilon ){
		step[ind] = 0;
		return;
	}// for a boundary node, do not change its value

	double p1,p2,p3,p4,p5,p6,p7;
	double r1,r2,r3,l1,l2,l3;
	double v_fore, v_back;

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

	v_fore = cpr[ind];
	v_back = cpl[ind];

	double xR, xL;
	weno_derivative_boundary(xR,xL,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dx,v_fore,v_back);

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

	v_fore = cpf[ind];
	v_back = cpb[ind];

	double yF, yB;
	weno_derivative_boundary(yF,yB,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dy,v_fore,v_back);

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

	v_fore = cpu[ind];
	v_back = cpd[ind];

	double zU, zD;
	weno_derivative_boundary(zU,zD,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dz,v_fore,v_back);
	
	step[ind] = (min2(0,vx[ind]) * xR + max2(0,vx[ind]) * xL + 
				 min2(0,vy[ind]) * yF + max2(0,vy[ind]) * yB + 
				 min2(0,vz[ind]) * zU + max2(0,vz[ind]) * zD ) * deltat[ind];

}








