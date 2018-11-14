/* shared_utilities.hpp */
#ifndef SHARED_UTILITIES_HPP
#define SHARED_UTILITIES_HPP

template<typename T> __device__ inline 
T max2(T x, T y) {return (x<y) ? y : x; }

template<typename T> __device__ inline 
T min2(T x, T y) {return (x<y) ? x : y; }

template<typename T> __device__ inline 
T min_mod(T x, T y) { return (x*y<0) ? 0.0 : (fabs(x)<fabs(y) ? x : y); }

template<typename T> __device__ inline
T min_abs(T x, T y) { return (fabs(x)<fabs(y)) ? x : y; }

template<typename T> __device__ inline
double sign(T x) { return (x>0) ? 1.0 : -1.0; }

// convert subindex to linear index
template<typename T> __device__ inline
T sub2ind(T row_idx, T col_idx, T pge_idx, T rows, T cols, T pges)
{
	T row_idxn = min2(rows-1, max2(0, row_idx));
	T col_idxn = min2(cols-1, max2(0, col_idx));
	T pge_idxn = min2(pges-1, max2(0, pge_idx));

	T ind = pge_idxn * rows * cols + col_idxn * rows + row_idxn;

	return ind;
}

/********************************************************************************/
/********************************************************************************/
// a Polynomial template that will accept any number of coefficients
// pn(x) = c0 + c1*x + c2*x^2 + ... cn*x^n
template<int n> class Poly
{
private:
	double _coef[n+1]; // c0,c1,...cn
public:
	__device__ inline Poly() {for(int i=0; i<=n; i++) _coef[i] = 0.0;}
	__device__ inline double operator()(double x);
	__device__ inline Poly<n> & setInterpolation(double const * x, double const * y);
};

template<int n> __device__ inline
double Poly<n>::operator()(double x)
{
	double val(0);
	for(int i=0; i<=n; i++) val += _coef[i]*pow(x,i);
	return val;
}

// set _coef from Lagrange's interpolation polynomial
template<int n> __device__ inline
Poly<n> & Poly<n>::setInterpolation(double const * x, double const * y)
{
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
	// c[i][:] is the _coefficients of the polynomial with zeros at xx[i][:]
	double c[n+1][n+1];
	for(int i=0; i<=n; i++){
		c[i][0] = 1.0;
		for(int j=1; j<=n; j++){
			c[i][j] = 0.0;
		}
	}
	for(int i=0; i<=n; i++){
		// _coefficient for c[i][:]
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
	// now collect results in _coef[:]
	for(int i=0; i<=n; i++){
		_coef[i] = 0.0;
		for(int j=0; j<=n; j++){
			_coef[i] += kk[j] * c[j][n-i];
		}
	}
	return *this;
}
/********************************************************************************/
/********************************************************************************/

// a template function that will find root of T(x) between x1 and x2 with a tolerance toler
// with Brent's methods
template<typename T> __device__ inline 
double zeroin(T func, double x1, double x2, double toler)
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

#endif
