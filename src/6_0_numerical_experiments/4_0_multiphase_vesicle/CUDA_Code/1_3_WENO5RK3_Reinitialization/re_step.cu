/*******************************************************************************
 * use weno derivative to calculate the numerical Hamiltonian for reinitialization
 * scheme
 * 5th order polynomial interpolation will be used to locate the boundary
 * weno53 scheme will be implemented on a nonuniform stencil near boundary 
 * central weno 6th order accurate scheme will be applied at nodes not immediately
 * next to the boundary
 ******************************************************************************/
#include "shared_utilities.hpp"

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
void select_stencil(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
		h0 = p4; x0 = 0.0;

		if(r1<ds){
			x1 = r1;
			x2 = ds;
			x3 = 2*ds;
			h1 = 0.0;
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
			h1m = 0.0;
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

__device__ inline
void weno_nonuniform(double & d_fore, double & d_back, double h3m, double h2m, double h1m, double h0, double h1, double h2, double h3, double x3m, double x2m, double x1m, double x0, double x1, double x2, double x3)
{
	// first divided differences, i.e. cell averages of derivatives
	double d1_2_5   = (h3  - h2)  / (x3  - x2) ;
	double d1_1_5   = (h2  - h1)  / (x2  - x1) ;
	double d1_0_5   = (h1  - h0)  / (x1  - x0) ;
	double d1_m0_5  = (h0  - h1m) / (x0  - x1m);
	double d1_m1_5  = (h1m - h2m) / (x1m - x2m);
	double d1_m2_5  = (h2m - h3m) / (x2m - x3m);

	// boundary nodes and cell averages
	double x_m2_5, x_m1_5, x_m0_5, x_0_5, x_1_5, x_2_5;
	double u_m2, u_m1, u_0, u_1, u_2;
	double v0, v1, v2, c0, c1, c2, s0, s1, s2;
	double epsilon = 1e-6, alpha0, alpha1, alpha2, sum, omega0, omega1, omega2;

	// calculate d_fore, choose a forward baised stencil
	x_m2_5 = x2m; x_m1_5 = x1m; x_m0_5 = x0; x_0_5 = x1; x_1_5 = x2; x_2_5 = x3;
	u_m2 = d1_m1_5; u_m1 = d1_m0_5; u_0 = d1_0_5; u_1 = d1_1_5; u_2 = d1_2_5;

	// now we calculate u_m0_5 from cell averages for different stencils
	v0 = u_1 
		+ (u_0 - u_1) * (1.0 + (x_0_5 - x_m0_5)/(x_1_5 - x_m0_5) + (x_0_5 - x_m0_5)/(x_2_5 - x_m0_5))
		+ (u_2 - u_1) * ((x_0_5 - x_m0_5)/(x_2_5 - x_m0_5)) * ((x_1_5 - x_m0_5)/(x_2_5 - x_0_5)) ;
	v1 = u_0 
		+ (u_m1 - u_0) * ((x_0_5 - x_m0_5)/(x_1_5 - x_m1_5)) * ((x_1_5 - x_m0_5)/(x_0_5 - x_m1_5))
		- (u_1  - u_0) * ((x_0_5 - x_m0_5)/(x_1_5 - x_m1_5)) * ((x_m0_5 - x_m1_5)/(x_1_5 - x_m0_5)) ;
	v2 = u_m1 
		+ (u_0  - u_m1) * ((x_m0_5 - x_m1_5)/(x_0_5 - x_m2_5)) * ((x_m0_5 - x_m2_5)/(x_0_5 - x_m1_5))
		- (u_m2 - u_m1) * ((x_m0_5 - x_m1_5)/(x_0_5 - x_m2_5)) * ((x_0_5 - x_m0_5)/(x_m0_5 - x_m2_5)) ;
	// optimal weights in smooth region
	c0 = ((x_m0_5 - x_m1_5)/(x_2_5 - x_m2_5)) * ((x_m0_5 - x_m2_5)/(x_2_5 - x_m1_5)) ;
	c1 = ((x_m0_5 - x_m2_5)/(x_2_5 - x_m2_5)) * ((x_2_5 - x_m0_5)/(x_2_5 - x_m1_5))
		* (1.0 + (x_2_5 - x_m1_5)/(x_1_5 - x_m2_5)) ;
	c2 = ((x_1_5 - x_m0_5)/(x_2_5 - x_m2_5)) * ((x_2_5 - x_m0_5)/(x_1_5 - x_m2_5)) ;
	// smoothness indicator
	{
	s0 = 4.0 * pow( (x_0_5 - x_m0_5)/(x_2_5 - x_m0_5), 2) * 
		(
			pow( (u_2 - u_1)/(x_2_5 - x_0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) + (x_1_5-x_m0_5) * (x_1_5 - x_0_5) 
				)
			+ (u_2 - u_1) * (u_0 - u_1) / ((x_2_5 - x_0_5)*(x_1_5 - x_m0_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5,2) + 2.0 * (x_1_5 - x_m0_5) * (x_1_5 - x_0_5) 
				  + (x_2_5 - x_m0_5) * (2.0 * x_1_5 - x_0_5 - x_m0_5) 
				)
			+ pow( (u_0 - u_1)/(x_1_5 - x_m0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) 
				  + (x_2_5 + x_1_5 - 2.0 * x_m0_5) * (x_2_5 + x_1_5 - x_0_5 - x_m0_5) 
				)
		);
	s1 = 4.0 * pow( (x_0_5 - x_m0_5)/(x_1_5 - x_m1_5), 2) * 
		(
			pow( (u_m1 - u_0)/(x_0_5 - x_m1_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) + (x_1_5 - x_m0_5) * (x_1_5 - x_0_5) 
				)
			+ (u_1 - u_0) * (u_m1 - u_0) / ((x_1_5 - x_m0_5)*(x_0_5 - x_m1_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5,2) - (x_1_5 - x_0_5)*(x_m0_5 - x_m1_5) 
				  - (x_1_5 - x_m0_5)*(x_0_5 - x_m1_5) 
				)
			+ pow( (u_1 -u_0)/(x_1_5 - x_m0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5, 2) + (x_m0_5 - x_m1_5)*(x_0_5 - x_m1_5)
				)
		);
	s2 = 4.0 *pow( (x_0_5 - x_m0_5)/(x_0_5 - x_m2_5), 2) * 
		(
			pow( (u_m2 - u_m1)/(x_m0_5 - x_m2_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5, 2) + (x_0_5 - x_m1_5)*(x_m0_5 - x_m1_5)
				)
			+ (u_0 - u_m1)*(u_m2 - u_m1) / ((x_0_5 - x_m1_5)*(x_m0_5 - x_m2_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5, 2)+ 2.0 * (x_0_5 - x_m1_5)*(x_m0_5 - x_m1_5) 
				  + (x_0_5 - x_m2_5)*(x_0_5 + x_m0_5 - 2.0 * x_m1_5)
				)
			+ pow( (u_0 - u_m1)/(x_0_5 - x_m1_5), 2) *
				( 10.0 * pow(x_0_5 - x_m0_5, 2) 
				  + (2.0 * x_0_5 - x_m2_5 - x_m1_5)*(x_0_5 + x_m0_5 - x_m1_5 - x_m2_5)
				)
		);
	}
	// optimal weights
	alpha0 = c0 / pow( (s0 + epsilon), 2);
	alpha1 = c1 / pow( (s1 + epsilon), 2);
	alpha2 = c2 / pow( (s2 + epsilon), 2);
	sum = alpha0 + alpha1 + alpha2;
	omega0 = alpha0 / sum;
	omega1 = alpha1 / sum;
	omega2 = alpha2 / sum;

	d_fore = v0 * omega0 + v1 * omega1 + v2 * omega2;

	// calculate d_back, choose a backward baised stencil
	x_m2_5 = x3m; x_m1_5 = x2m; x_m0_5 = x1m; x_0_5 = x0; x_1_5 = x1; x_2_5 = x2;
	u_m2 = d1_m2_5; u_m1 = d1_m1_5; u_0 = d1_m0_5; u_1 = d1_0_5; u_2 = d1_1_5; 

	// now we calculate u_0_5 from cell averages for different stencils
	v0 = u_1 
		+ (u_0 - u_1) * ((x_1_5 - x_0_5)/(x_2_5 - x_m0_5)) * ((x_2_5 - x_0_5)/(x_1_5 - x_m0_5))
		- (u_2 - u_1) * ((x_1_5 - x_0_5)/(x_2_5 - x_m0_5)) * ((x_0_5 - x_m0_5)/(x_2_5 - x_0_5)) ;
	v1 = u_0
		+ (u_1 - u_0)  * ((x_0_5 - x_m0_5)/(x_1_5 - x_m1_5)) * ((x_0_5 - x_m1_5)/(x_1_5 - x_m0_5))
		- (u_m1 - u_0) * ((x_0_5 - x_m0_5)/(x_1_5 - x_m1_5)) * ((x_1_5 - x_0_5)/(x_0_5 - x_m1_5)) ;
	v2 = u_m1
		+ (u_m2 - u_m1) * ((x_0_5 - x_m0_5)/(x_m0_5 - x_m2_5)) * ((x_0_5 - x_m1_5)/(x_0_5 - x_m2_5))
		+ (u_0 - u_m1) * (1.0 + (x_0_5 - x_m0_5)/(x_0_5 - x_m1_5) + (x_0_5 - x_m0_5)/(x_0_5 - x_m2_5)) ;
	// optimal weights in smooth region
	c0 = ((x_0_5 - x_m2_5)/(x_2_5 - x_m2_5)) * ((x_0_5 - x_m1_5)/(x_2_5 -x_m1_5)) ;
	c1 = ((x_0_5 - x_m2_5)/(x_2_5 - x_m2_5)) * ((x_2_5 - x_0_5)/(x_2_5 - x_m1_5))
		* (1.0 + (x_2_5 - x_m1_5)/(x_1_5 - x_m2_5)) ;
	c2 = ((x_1_5 - x_0_5)/(x_2_5 - x_m2_5)) * ((x_2_5 - x_0_5)/(x_1_5 - x_m2_5)) ;
	// smoothness indicator
	{
	s0 = 4.0 * pow( (x_0_5 - x_m0_5)/(x_2_5 - x_m0_5), 2) * 
		(
			pow( (u_2 - u_1)/(x_2_5 - x_0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) + (x_1_5-x_m0_5) * (x_1_5 - x_0_5) 
				)
			+ (u_2 - u_1) * (u_0 - u_1) / ((x_2_5 - x_0_5)*(x_1_5 - x_m0_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5,2) + 2.0 * (x_1_5 - x_m0_5) * (x_1_5 - x_0_5) 
				  + (x_2_5 - x_m0_5) * (2.0 * x_1_5 - x_0_5 - x_m0_5) 
				)
			+ pow( (u_0 - u_1)/(x_1_5 - x_m0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) 
				  + (x_2_5 + x_1_5 - 2.0 * x_m0_5) * (x_2_5 + x_1_5 - x_0_5 - x_m0_5) 
				)
		);
	s1 = 4.0 * pow( (x_0_5 - x_m0_5)/(x_1_5 - x_m1_5), 2) * 
		(
			pow( (u_m1 - u_0)/(x_0_5 - x_m1_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5,2) + (x_1_5 - x_m0_5) * (x_1_5 - x_0_5) 
				)
			+ (u_1 - u_0) * (u_m1 - u_0) / ((x_1_5 - x_m0_5)*(x_0_5 - x_m1_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5,2) - (x_1_5 - x_0_5)*(x_m0_5 - x_m1_5) 
				  - (x_1_5 - x_m0_5)*(x_0_5 - x_m1_5) 
				)
			+ pow( (u_1 -u_0)/(x_1_5 - x_m0_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5, 2) + (x_m0_5 - x_m1_5)*(x_0_5 - x_m1_5)
				)
		);
	s2 = 4.0 *pow( (x_0_5 - x_m0_5)/(x_0_5 - x_m2_5), 2) * 
		(
			pow( (u_m2 - u_m1)/(x_m0_5 - x_m2_5), 2) * 
				( 10.0 * pow(x_0_5 - x_m0_5, 2) + (x_0_5 - x_m1_5)*(x_m0_5 - x_m1_5)
				)
			+ (u_0 - u_m1)*(u_m2 - u_m1) / ((x_0_5 - x_m1_5)*(x_m0_5 - x_m2_5)) *
				( 20.0 * pow(x_0_5 - x_m0_5, 2)+ 2.0 * (x_0_5 - x_m1_5)*(x_m0_5 - x_m1_5) 
				  + (x_0_5 - x_m2_5)*(x_0_5 + x_m0_5 - 2.0 * x_m1_5)
				)
			+ pow( (u_0 - u_m1)/(x_0_5 - x_m1_5), 2) *
				( 10.0 * pow(x_0_5 - x_m0_5, 2) 
				  + (2.0 * x_0_5 - x_m2_5 - x_m1_5)*(x_0_5 + x_m0_5 - x_m1_5 - x_m2_5)
				)
		);
	}
	// optimal weights
	alpha0 = c0 / pow( (s0 + epsilon), 2);
	alpha1 = c1 / pow( (s1 + epsilon), 2);
	alpha2 = c2 / pow( (s2 + epsilon), 2);
	sum = alpha0 + alpha1 + alpha2;
	omega0 = alpha0 / sum;
	omega1 = alpha1 / sum;
	omega2 = alpha2 / sum;

	d_back = v0 * omega0 + v1 * omega1 + v2 * omega2;
}

// calculate weno derivative at p4: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// where px are level set function values at node x
// lx, rx are distance to the left/right node 
__device__ inline
void weno_derivative_boundary(double & d_fore, double & d_back, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
	bool cross_interface = p3*p4<0 || p4*p5<0;

	if(cross_interface){
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		weno_nonuniform(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3); 
	}// for nodes IMMEDIATELY adjacent to the boundary, use cubic ENO interpolant
	else{
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p5) / ds;
		double v6 = (p7 - p6) / ds;
		d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
	}// if not a node IMMEDIATELY adjacent to the boundary, calculate weno derivatives as usual
}

__global__
void re_step(double * step, double const * lsf, bool const * mask, double const * deltat, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	double epsilon = 1e-6 * dx;
	if( xpr[ind]< epsilon || xpl[ind]<epsilon || ypf[ind]<epsilon || ypb[ind]<epsilon || zpu[ind]<epsilon || zpd[ind]<epsilon ){
		step[ind] = 0;
		return;
	}// for a boundary node, do not change its value

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

	double xR, xL;
	weno_derivative_boundary(xR,xL,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dx);

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

	double yF, yB;
	weno_derivative_boundary(yF,yB,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dy);

	int upup1	= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	int upup2	= sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
	int upup3	= sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);
	int down1	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2	= sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3	= sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

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

	double zU, zD;
	weno_derivative_boundary(zU,zD,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dz);

	if (mask[ind]) {
		step[ind] = ( sqrt(	max2(pow(min2(0.0,xL),2),pow(max2(0.0,xR),2)) + 
							max2(pow(min2(0.0,yB),2),pow(max2(0.0,yF),2)) + 
							max2(pow(min2(0.0,zD),2),pow(max2(0.0,zU),2)) ) - 1
					) * deltat[ind] * (-1.);
	} else{
		step[ind] = ( sqrt(	max2(pow(max2(0.0,xL),2),pow(min2(0.0,xR),2)) + 
							max2(pow(max2(0.0,yB),2),pow(min2(0.0,yF),2)) + 
							max2(pow(max2(0.0,zD),2),pow(min2(0.0,zU),2)) ) - 1
					) * deltat[ind] * (1.);
	}
}

