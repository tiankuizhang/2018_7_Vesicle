/*******************************************************************************
 * make corrections to xpr etc.
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
	return (x*y<0) ? 0 : (fabs(x)<fabs(y) ? x : y);
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
	//int row_idxn = (row_idx + rows) % rows;
	//int col_idxn = (col_idx + cols) % cols;
	//int pge_idxn = (pge_idx + pges) % pges;

	int row_idxn = min2(rows-1, max2(0, row_idx));
	int col_idxn = min2(cols-1, max2(0, col_idx));
	int pge_idxn = min2(pges-1, max2(0, pge_idx));


	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

/****************************************************************************** 
 * calculate distance to the bundary. 
 * if v0*f2<0, return distance from v0 to boundary
 ******************************************************************************/
__device__ inline
double sp(double v1, double v0, double v2, double v3, double ds)
{
	double epsilon=10e-10;
	double p2l = v2 - 2.0 * v0 + v1;
	double p2r = v3 - 2.0 * v2 + v0;
	double p2 = min_mod(p2l, p2r);
	if(p2>epsilon || p2<-epsilon){
		double disc = pow((0.5*p2-v0-v2),2) - 4.*v0*v2;
		double dist =  ds * (0.5 + (v0-v2 - sign(v0-v2)*sqrt(disc)) / p2 );
		return dist;
	}else{
		return ( ds * v0 / (v0 - v2) ); 
	}

}

// given 4 points (0,v0),(s,v1),(2*s,v2),(3*s,v3) and v1*v2<0
// calculate coefficients of the cubic interpolant c0,c1,c2,c3
// p3(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3;
__device__ inline
double cubic_distance(double v0, double v1, double v2, double v3, double s)
{
	// calculate the interpolant coefficient
	double c0 = v0;
 	double c1 = (     3.0 * (v1-v0) - 3.0/2.0 * (v2-v0) + 1.0/3.0 * (v3-v0) ) / s;
 	double c2 = (-5.0/2.0 * (v1-v0) +     2.0 * (v2-v0) - 1.0/2.0 * (v3-v0) ) / pow(s,2);
 	double c3 = ( 1.0/2.0 * (v1-v0) - 1.0/2.0 * (v2-v0) + 1.0/6.0 * (v3-v0) ) / pow(s,3);
	/* It is EXTREMELY important to use float point numbers 1.0/2.0 instead of 1/2 
	 * the latter will give (double)(int)(1/2) = 0.0 instead of 0.5
	 */

	// now use Newton's method to find root
	double xc = s + s * v1 / (v1 - v2); // initial guess

	int iter = 0;
	int const max_iter = 50;
	double const max_error = 1e-14;
	double diff = 1;

	while( diff>max_error && iter<max_iter){
		// Newton's method
		double f = c0 + c1 * xc + c2 * xc*xc + c3 * xc*xc*xc;
		double d = c1 + 2 * c2 * xc + 3 * c3 * xc * xc;
		double new_xc = xc - f / d;

		diff = fabs (f / d);
		iter++;

		xc = new_xc;
	}

	return (xc - s);
}

/****************************************************************************** 
 * calculate Eno derivatives at node v0: [v4,v1,v0,v2,v3]
 ******************************************************************************/
__device__ inline
void eno_derivative(double & sR, double & sL, double v4, double v1, double v0, double v2, double v3, double pr, double pl, double ds)
{
	double p2m;

	double p2 = v1 - 2.0 * v0 + v2;

	double p2r = v0 - 2.0 * v2 + v3;
	p2m = 0.5 * min_mod(p2, p2r) / pow(ds, 2);
	double vr = (pr==ds) ? v2 : 0;
	sR = (vr - v0) / pr - pr * p2m;

	double p2l = v0 - 2.0 * v1 + v4;
	p2m = 0.5 * min_mod(p2, p2l) / pow(ds, 2);
	double vl = (pl==ds) ? v1 : 0;
	sL = (v0 - vl) / pl + pl * p2m;

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
void select_stencil_backup(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
		h0 = p4; x0 = 0.0;

		double hr[6] = {0,0,0,0,0,0};
		double xr[6] = {0,0,0,0,0,0};
		int i = 0;

		if(r1<ds){
			xr[i] = r1; hr[i] = 0.0; i++;
		}else{
			xr[i] = ds;	hr[i] = p5; i++;
		}
		if(r2<ds){
			xr[i] = r2+ds; hr[i] = 0.0; i++;
		}else{
			xr[i] = 2*ds; hr[i] = p6; i++;
		}
		if(r3<ds){
			xr[i] = r3+2*ds; hr[i] = 0.0; i++;
		}else{
			xr[i] = 3*ds; hr[i] = p7; i++;
		}

		x1 = xr[0]; x2 = xr[1]; x3 = xr[2];
		h1 = hr[0]; h2 = hr[1]; h3 = hr[2];

		double hl[6] = {0,0,0,0,0,0};
		double xl[6] = {0,0,0,0,0,0};
		i = 0;

		if(l1<ds){
			xl[i] = l1; hl[i] = 0.0; i++;
		}else{
			xl[i] = ds;	hl[i] = p3; i++;
		}
		if(l2<ds){
			xl[i] = l2+ds; hl[i] = 0.0; i++;
		}else{
			xl[i] = 2*ds; hl[i] = p2; i++;
		}
		if(l3<ds){
			xl[i] = l3+2*ds; hl[i] = 0.0; i++;
		}else{
			xl[i] = 3*ds; hl[i] = p1; i++;
		}

		h1m =  hl[0]; h2m =  hl[1]; h3m =  hl[2];
		x1m = -xl[0]; x2m = -xl[1]; x3m = -xl[2];
}

// given a stencil across the boundary: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// create a new stencil (x3m,h3m),(x2m,h2m),(x1m,h1m),(x0,h0),(x1,h1),(x2,h2),(x3,h3) including boundary nodes
__device__ inline
void select_stencil(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
		h0 = p4; x0 = 0.0;

		double hr[6] = {0,0,0,0,0,0};
		double xr[6] = {0,0,0,0,0,0};
		int i = 0;

		if(r1<ds){
			xr[i] = r1; hr[i] = 0.0; i++;
			xr[i] = ds;	hr[i] = p5; i++;
		}else{
			xr[i] = ds;	hr[i] = p5; i++;
		}
		if(r2<ds){
			xr[i] = r2+ds; hr[i] = 0.0; i++;
			xr[i] = 2*ds; hr[i] = p6; i++;
		}else{
			xr[i] = 2*ds; hr[i] = p6; i++;
		}
		if(r3<ds){
			xr[i] = r3+2*ds; hr[i] = 0.0; i++;
			xr[i] = 3*ds; hr[i] = p7; i++;
		}else{
			xr[i] = 3*ds; hr[i] = p7; i++;
		}

		x1 = xr[0]; x2 = xr[1]; x3 = xr[2];
		h1 = hr[0]; h2 = hr[1]; h3 = hr[2];

		double hl[6] = {0,0,0,0,0,0};
		double xl[6] = {0,0,0,0,0,0};
		i = 0;

		if(l1<ds){
			xl[i] = l1; hl[i] = 0.0; i++;
			xl[i] = ds;	hl[i] = p3; i++;
		}else{
			xl[i] = ds;	hl[i] = p3; i++;
		}
		if(l2<ds){
			xl[i] = l2+ds; hl[i] = 0.0; i++;
			xl[i] = 2*ds; hl[i] = p2; i++;
		}else{
			xl[i] = 2*ds; hl[i] = p2; i++;
		}
		if(l3<ds){
			xl[i] = l3+2*ds; hl[i] = 0.0; i++;
			xl[i] = 3*ds; hl[i] = p1; i++;
		}else{
			xl[i] = 3*ds; hl[i] = p1; i++;
		}

		h1m =  hl[0]; h2m =  hl[1]; h3m =  hl[2];
		x1m = -xl[0]; x2m = -xl[1]; x3m = -xl[2];
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
void weno_derivative_boundary_backup(double & d_fore, double & d_back, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
	bool out_side = p1>0 && p2>0 && p3>0 && p4>0 && p5>0 && p6>0 && p7>0;
	bool in_side  = p1<0 && p2<0 && p3<0 && p4<0 && p5<0 && p6<0 && p7<0;

	if(out_side || in_side){
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p4) / ds;
		double v6 = (p7 - p6) / ds;
		d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
	}// if not a boundary node, calculate weno derivatives as usual
	else{
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
	}
	
}

// calculate weno derivative at p4: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// where px are level set function values at node x
// lx, rx are distance to the left/right node 
__device__ inline
void weno_derivative_boundary(double & d_fore, double & d_back, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{

	bool cross_interface = p3*p4<0 || p4*p5<0;

	if(!cross_interface){
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p5) / ds;
		double v6 = (p7 - p6) / ds;
		d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
	}// if not a boundary node, calculate weno derivatives as usual
	else{
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
	}
	
}

// make corrections to xpr etc
__global__
void boundary_correction(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	double f2;
	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);
	double f0 = lsf[ind];

	int right = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	f2 = lsf[right];
	if(f0*f2<0){
		int left = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
		int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
		xpr[ind] = cubic_distance(lsf[left], f0, f2, lsf[right2], dx);
		xpl[right] = dx - xpr[ind];
	}

	int front = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	f2 = lsf[front];
	if(f0*f2<0){
		int back = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
		int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
		ypf[ind] = cubic_distance(lsf[back], f0, f2, lsf[front2], dy);
		ypb[front] = dy - ypf[ind];
	}

	int up = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	f2 = lsf[up];
	if(f0*f2<0){
		int down = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
		int up2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
		zpu[ind] = cubic_distance(lsf[down], f0, f2, lsf[up2], dz);
		zpd[up] = dz - zpu[ind];
	}

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
	//eno_derivative(xR,xL,p2,p3,p4,p5,p6,r1,l1,dx);
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
	//eno_derivative(yF,yB,p2,p3,p4,p5,p6,r1,l1,dy);
	weno_derivative_boundary(yF,yB,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dy);

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

	double zU, zD;
	//eno_derivative(zU,zD,p2,p3,p4,p5,p6,r1,l1,dz);
	weno_derivative_boundary(zU,zD,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,dz);

	
	if (mask[ind]) {
		step[ind] = ( sqrt(	max2(pow(min2(0,xL),2),pow(max2(0,xR),2)) + 
							max2(pow(min2(0,yB),2),pow(max2(0,yF),2)) + 

							max2(pow(min2(0,zD),2),pow(max2(0,zU),2)) ) - 1)
					* deltat[ind] * (-1.);
	} else{
		step[ind] = ( sqrt(	max2(pow(max2(0,xL),2),pow(min2(0,xR),2)) + 
							max2(pow(max2(0,yB),2),pow(min2(0,yF),2)) + 
							max2(pow(max2(0,zD),2),pow(min2(0,zU),2)) ) - 1)
			* deltat[ind] * (1.);
	}
}

