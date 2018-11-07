/*******************************************************************************
 * use weno derivative to calculate the numerical Hamiltonian for reinitialization
 * scheme
 * 5th order polynomial interpolation will be used to locate the boundary
 * weno53 scheme will be implemented on a nonuniform stencil near boundary 
 * central weno 6th order accurate scheme will be applied at nodes not immediately
 * next to the boundary
 ******************************************************************************/

namespace su // shared utilities
{
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
	int row_idxn = min2(rows-1, max2(0, row_idx));
	int col_idxn = min2(cols-1, max2(0, col_idx));
	int pge_idxn = min2(pges-1, max2(0, pge_idx));

	int ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}
}

namespace bd // boundary distance
{
	using namespace su;
// a Polynomial template that will accept any number of coefficients
// pn(x) = c0 + c1*x + c2*x^2 + ... cn*x^n
template<int n> class Poly
{
private:
	double coef[n+1]; // c0,c1,...cn
public:
	__device__ inline
	Poly<n>(void) {
		for(int i=0; i<=n; i++) coef[i] = 0.0;
	}
	__device__ inline
	int setCoefficient(double const * c){
		for(int i=0; i<=n; i++) coef[i] = c[i];
		return 0;
	}
	// set coef from Lagrange's interpolation polynomial
	__device__ inline
	int setInterpolation(double const * x, double const * y){
		// n+1 points are needed to determine a nth order poly
		double kk[n+1]; // multiplying factor for each Lagrange interpolation term
		double xx[n+1][n]; // corresponding zeros of each Lagrange interpolation term
		for(int i=0; i<=n; i++){
			kk[i] = y[i];
			int ind(0);
			for(int j=0; j<=n; j++){
				if(j!=i){
					kk[i] /= x[i]-x[j];
					xx[i][ind] = x[j];
					ind++;
				}
			}
		}
		// c[i][:] is the coefficients of the polynomial with zeros at xx[i][:]
		double c[n+1][n+1];
		for(int i=0; i<=n; i++){
			c[i][0] = 1.0;
			for(int j=1; j<=n; j++){
				c[i][j] = 0.0;
			}
		}
		for(int i=0; i<=n; i++){
			// coefficient for c[i][:]
			for(int j=0; j<n; j++){
				// n iteration for n points
				double tmp[n+1];
				for(int k=1; k<=j+1; k++){
					// expand recursion formula
					tmp[k] = c[i][k] - xx[i][j]*c[i][k-1];
				}
				for(int k=1; k<=j+1; k++){
					c[i][k] = tmp[k];
				}
			}
		}
		// now collect results in coef[:]
		for(int i=0; i<=n; i++){
			coef[i] = 0.0;
			for(int j=0; j<=n; j++){
				coef[i] += kk[j] * c[j][n-i];
			}
		}
		return 0;
	}
	__device__ inline
	double operator()(double x){
		double val(0);
		for(int i=0; i<=n; i++) val += coef[i]*pow(x,i);
		return val;
	}
};

// a template function that will find root of T(x) between x1 and x2 with a tolerance toler
// with Brent's methods
template<typename T>
__device__ inline double zeroin(T func, double x1, double x2, double toler)
{
	int const ITMAX = 100; // Maximum allowed number of iterations

	double d,e,p,q,r,s,xm;
	double a=x1, b=x2, c=x2;
	double fa=func(a), fb=func(b), fc=fb;

	for(int iter=0; iter<ITMAX; iter++){
		// make sure result be bound by b,c; b is the best result; a is the previous result
		if((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)){
			c = a; fc = fa;
			d = b - a; e = d;
		}
		if(fabs(fc) < fabs(fb)){
			a = b; b = c; c = a;
			fa = fb; fb = fc; fc = fa;
		}
		// convergence check
		xm = 0.5*(c-b);
		if(fabs(xm) <= toler || fb == 0.0) return b;
		// choose bisection or interpolation
		if(fabs(e) < toler || abs(fa) <= abs(fb)){
			// bounds decreasing too slowly, use bisection
			d = xm; e = xm;
		}else{
			// attempt interpolation
			s = fb/fa;
			if(a == c){
				// linear interpolation
				p = 2.0 * xm * s;
				q = 1.0 - s;
			}else{
				// inverse quadratice interpolation
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			}
			if (p > 0) q = -q;
			p = fabs(p);
			// is interpolation acceptable
			double Min1 = 3.0 * xm * q - abs(toler*q);
			double Min2 = abs(e*q);
			if(2.0*p < min2(Min1,Min2)){
				// accept interpolation
				e = d; d = p/q;
			}else{
				// use bisection
				d = xm; e = d;
			}
		}
		// move last best guess to a
		a = b; fa = fb;
		// elvaluate new trial root
		if(fabs(d) > toler){
			b = b + d;
		}else{
			b += sign(xm) * toler;
		}
		fb = func(b);
	}
	return b;
};

/****************************************************************************** 
 * calculate distance to the bundary with qudratic eno 
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
	// IMPORTANT
	bool kink = v0*v1 <0 || v2*v3 < 0;
	if(kink){
		return sp(v0,v1,v2,v3,s);
	}// if near a kink, use cubic ENO interpolant

	// create a cubic interpolation
	double xx[4] = {0.0,s,2*s,3*s};
	double yy[4] = {v0,v1,v2,v3};
	Poly<3> p3;
	p3.setInterpolation(xx,yy);

	// find zero of the polynomial
	double xc = zeroin<Poly<3> >(p3,s,2*s,3e-16);

	// result should lie between s and 2*s
	xc = max2(xc,s);
	xc = min2(xc,2*s);

	return (xc - s);
}

// calculate distance to boundary with a 6th order accurate approximation
// i.e., with a 5th order polynomial
// given 6 points (0,v0),(s,v1),(2s,v2),(3s,v3),(4s,v4),(5s,v5) and v2*v3<0
__device__ inline
double six_distance(double v0, double v1, double v2, double v3, double v4, double v5, double s)
{
	// IMPORTANT
	// if near a kink, use cubic ENO interpolant
	bool kink1 = v1*v2 < 0 || v3*v4 < 0;
	if(kink1) return sp(v1,v2,v3,v4,s);

	// if one grid away from a kink use cubic interpolation
	//bool kink2 = v0*v1 < 0 || v4*v5 < 0;
	//if(kink2 && false) return cubic_distance(v1,v2,v3,v3,s);
	// should not be used due to unexpected degradation of accuracy

	// create a cubic interpolation
	double xx[6] = {0.0,s,2*s,3*s,4*s,5*s};
	double yy[6] = {v0,v1,v2,v3,v4,v5};
	Poly<5> p5;
	p5.setInterpolation(xx,yy);

	// find zero of the polynomial
	double xc = zeroin<Poly<5> >(p5,2*s,3*s,3e-16);

	// result should lie between 2*s and 3*s
	xc = max2(xc,2*s);
	xc = min2(xc,3*s);

	return (xc - 2*s);
}
}

namespace rs // reinitialization step
{
	using namespace su;
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

__device__ inline
double weno6_onesided_derivative(double vm2, double vm1, double v0, double v1, double v2, double v3)
{
	// different choices of ENO derivatives
	double phi0 = (2.0  * vm2 - 7.0 * vm1 + 11.0 * v0) / 6.0;
	double phi1 = (      -vm1 + 5.0 * v0  + 2.0  * v1) / 6.0;
	double phi2 = (2.0  * v0  + 5.0 * v1  -        v2) / 6.0;
	double phi3 = (11.0 * v1  - 7.0 * v2  + 2.0  * v3) / 6.0;

	// smoothness indicator
	double S0 = 1./4. * pow(vm2 - 4.*vm1 + 3.*v0, 2) + 13./12. * pow(vm2 - 2.*vm1 + v0, 2);
	double S1 = 1./4. * pow(vm1 - v1, 2)             + 13./12. * pow(vm1 - 2.*v0  + v1, 2);
	double S2 = 1./4. * pow(3.*v0 - 4*v1 + v2, 2)    + 13./12. * pow(v0  - 2.*v1  + v2, 2);
	double S6 = 
		(
			271779. * vm2 * vm2
			+ vm2 * ( 2380800.  * vm1 + 4086352.  * v0 - 3462252.  * v1 + 1458762. * v2 - 245620.  * v3)
			+ vm1 * ( 5653317.  * vm1 - 20427884. * v0 + 17905032. * v1 - 7727988. * v2 + 1325006. * v3)
			+ v0  * ( 19510972. * v0  - 35817664. * v1 + 15929912. * v2 - 2792660. * v3)
			+ v1  * ( 17195652. * v1  - 15880404. * v2 + 2863984.  * v3)
			+ v2  * ( 3824847.  * v2  - 1429976.  * v3)
			+ 139633. * v3 * v3
		) / 10080.;
	double S3 = S6;
	double tau6 = S6 - (S0 + S2 + 4.0 * S1) / 6.0;

	double epsilon = 1e-40;
	//double epsilon = 1e-6;
	double C = 20.0;
	double d0 = 1./20., d1 = 9./20., d2 = 9./20., d3 = 1./20.;

	double alpha0 = d0 * (C + tau6 / (S0 + epsilon) );
	double alpha1 = d1 * (C + tau6 / (S1 + epsilon) );
	double alpha2 = d2 * (C + tau6 / (S2 + epsilon) );
	double alpha3 = d3 * (C + tau6 / (S3 + epsilon) );

	double sum = alpha0 + alpha1 + alpha2 + alpha3;
	double omega0 = alpha0 / sum;
	double omega1 = alpha1 / sum;
	double omega2 = alpha2 / sum;
	double omega3 = alpha3 / sum;

	return (omega0 * phi0 + omega1 * phi1 + omega2 * phi2 + omega3 * phi3);
}

// given a stencil across the boundary: p1<-l3-p2<-l2-p3<-l1-p4-r1->p5-r2->p6-r3->p7
// create a new stencil (x3m,h3m),(x2m,h2m),(x1m,h1m),(x0,h0),(x1,h1),(x2,h2),(x3,h3) including boundary nodes
__device__ inline
void select_stencil_1(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
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
void select_stencil_2(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
		h0 = p4; x0 = 0.0;

		if(r2<ds){
			x1 = ds;
			x2 = ds + r2;
			x3 = 2*ds;
			h1 = p5;
			h2 = 0.0;
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
			x1m = - ds;
			x2m = - ds - l1;
			x3m = - 2*ds;
			h1m = p3;
			h2m = 0.0;
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
void select_stencil_0(double & h3m, double & h2m, double & h1m, double & h0, double & h1, double & h2, double & h3, double & x3m, double & x2m, double & x1m, double & x0, double & x1, double & x2, double & x3, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
		h0 = p4; x0 = 0.0;
		double r[3], p[3], x[7], h[7];
		int i, j;
		x[0] = x0; h[0] = h0;
		// set (x1,h1),(x2,h2),(x3,h3) from r1,r2,r3 and p5, p6, p7
		r[0] = r1; r[1] = r2; r[2] = r3; p[0] = p5; p[1] = p6; p[2] = p3;
		i = 0, j = 1;
		while(i < 3){
			if(r[i] < ds){
				x[j] = x[j-1] + r[i]; x[j+1] = x[j-1] + ds; 
				h[j] = 0.0; h[j+1] = p[i];
				j += 2; i++;
			}else{
				x[j] = x[j-1] + ds; h[j] = p[i]; j++; i++;
			}
		}
		x1 = x[1]; x2 = x[2]; x3 = x[3];
		h1 = h[1]; h2 = h[2]; h3 = h[3];
		// set (x1m,h1m),(x2m,h2m),(x3m,h3m) from l1,l2,l3 and p3,p2,p1
		r[0] = l1; r[1] = l2; r[2] = l3; p[0] = p3; p[1] = p2; p[2] = p1;
		i = 0; j = 1;
		while(i < 3){
			if(r[i] < ds){
				x[j] = x[j-1] - r[i]; x[j+1] = x[j-1] - ds;
				h[j] = 0.0; h[j+1] = p[i];
				j += 2; i++;
			}else{
				x[j] = x[j-1] - ds; h[j] = p[i]; j++; i++;
			}
		}
		x1m = x[1]; x2m = x[2]; x3m = x[3];
		h1m = h[1]; h2m = h[2]; h3m = h[3];



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

__device__ inline
void weno_nonuniform(double & d_fore, double & d_back, double h3m, double h2m, double h1m, double h0, double h1, double h2, double h3, double x3m, double x2m, double x1m, double x0, double x1, double x2, double x3)
{
	// first divided differences, i.e. cell averages
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
				( 20.0 * pow(x_0_5-x_m0_5, 2)+ 2.0 * (x_0_5 - x_1_5)*(x_m0_5 - x_m1_5) 
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
		+ (u_1 - u_0) * ((x_0_5 - x_m0_5)/(x_1_5 - x_m1_5)) * ((x_0_5 - x_m1_5)/(x_1_5 - x_m0_5))
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
				( 20.0 * pow(x_0_5-x_m0_5, 2)+ 2.0 * (x_0_5 - x_1_5)*(x_m0_5 - x_m1_5) 
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
	//bool cross_interface = p3*p4<0 || p4*p5<0 || p2*p3<0 || p5*p6<0;

	if(cross_interface){
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil_1(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		//select_stencil_0(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		//ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
		weno_nonuniform(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3); 
		// weno_nonuniform gives much better results
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
		//d_back = weno6_onesided_derivative(v1,v2,v3,v4,v5,v6); // less accurate
		//d_fore = weno6_onesided_derivative(v6,v5,v4,v3,v2,v1);
	}// if not a node IMMEDIATELY adjacent to the boundary, calculate weno derivatives as usual
}

__device__ inline
void weno_derivative_boundary_backup(double & d_fore, double & d_back, double p1, double p2, double p3, double p4, double p5, double p6, double p7, double r1, double r2, double r3, double l1, double l2, double l3, double ds)
{
	bool cross_interface_1 = p3*p4<0 || p4*p5<0;
	bool cross_interface_2 = p2*p3<0 || p5*p6<0;
	//bool cross_interface = p3*p4<0 || p4*p5<0 || p2*p3<0 || p5*p6<0;

	if(cross_interface_1){
		double h3m,h2m,h1m,h0,h1,h2,h3;
		double x3m,x2m,x1m,x0,x1,x2,x3;
		select_stencil_1(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		//ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
		weno_nonuniform(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
	}// for nodes IMMEDIATELY adjacent to the boundary, use cubic ENO interpolant
	else if(cross_interface_2){
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p5) / ds;
		double v6 = (p7 - p6) / ds;
		d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
		//d_back = weno6_onesided_derivative(v1,v2,v3,v4,v5,v6);
		//d_fore = weno6_onesided_derivative(v6,v5,v4,v3,v2,v1);
		//double h3m,h2m,h1m,h0,h1,h2,h3;
		//double x3m,x2m,x1m,x0,x1,x2,x3;
		//select_stencil_2(h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3,p1,p2,p3,p4,p5,p6,p7,r1,r2,r3,l1,l2,l3,ds);
		//ENO_cubic_derivative(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
		//weno_nonuniform(d_fore,d_back,h3m,h2m,h1m,h0,h1,h2,h3,x3m,x2m,x1m,x0,x1,x2,x3);
	}
	else{
		double v1 = (p2 - p1) / ds;
		double v2 = (p3 - p2) / ds;
		double v3 = (p4 - p3) / ds;
		double v4 = (p5 - p4) / ds;
		double v5 = (p6 - p5) / ds;
		double v6 = (p7 - p6) / ds;
		//d_back = weno_onesided_derivative(v1,v2,v3,v4,v5);
		//d_fore = weno_onesided_derivative(v6,v5,v4,v3,v2);
		d_back = weno6_onesided_derivative(v1,v2,v3,v4,v5,v6);
		d_fore = weno6_onesided_derivative(v6,v5,v4,v3,v2,v1);
	}// if not a node IMMEDIATELY adjacent to the boundary, calculate weno derivatives as usual
}
}

// make corrections to xpr etc with cubic interpolation
__global__
void boundary_correction_backup(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
{
	using namespace bd;

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

// make corrections to xpr etc with 5th order interpolation
// seems that 6th order interpolation gives the same results as cubic interpolation
__global__
void boundary_correction(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
{
	using namespace bd;

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
		int left2 = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
		int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
		int right3 = sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
		xpr[ind] = six_distance(lsf[left2], lsf[left], f0, f2, lsf[right2], lsf[right3], dx);
		xpl[right] = dx - xpr[ind];
	}

	int front = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	f2 = lsf[front];
	if(f0*f2<0){
		int back = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
		int back2 = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
		int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
		int front3 = sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);
		ypf[ind] = six_distance(lsf[back2], lsf[back], f0, f2, lsf[front2], lsf[front3], dy);
		ypb[front] = dy - ypf[ind];
	}

	int up = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	f2 = lsf[up];
	if(f0*f2<0){
		int down = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
		int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
		int up2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
		int up3 = sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);
		zpu[ind] = six_distance(lsf[down2], lsf[down], f0, f2, lsf[up2], lsf[up3], dz);
		zpd[up] = dz - zpu[ind];
	}

}

__global__
void re_step(double * step, double const * lsf, bool const * mask, double const * deltat, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	using namespace rs;

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

	//double smoothSign = p4 / sqrt(p4*p4 + sqrt(dx) );
	//double smoothSign = p4 / sqrt(p4*p4 + dx * dx );
	double smoothSign = sign(p4); // performs better than the versions above
	
	if (mask[ind]) {
		step[ind] = ( sqrt(	max2(pow(min2(0,xL),2),pow(max2(0,xR),2)) + 
							max2(pow(min2(0,yB),2),pow(max2(0,yF),2)) + 
							max2(pow(min2(0,zD),2),pow(max2(0,zU),2)) ) - 1
					) * deltat[ind] * smoothSign;
					//) * deltat[ind] * (-1.);
	} else{
		step[ind] = ( sqrt(	max2(pow(max2(0,xL),2),pow(min2(0,xR),2)) + 
							max2(pow(max2(0,yB),2),pow(min2(0,yF),2)) + 
							max2(pow(max2(0,zD),2),pow(min2(0,zU),2)) ) - 1
					) * deltat[ind] * smoothSign;
					//) * deltat[ind] * (1.);
	}
}

