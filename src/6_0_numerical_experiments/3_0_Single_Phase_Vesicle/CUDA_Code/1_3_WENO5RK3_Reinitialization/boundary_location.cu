#include "shared_utilities.hpp"
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

// make corrections to xpr etc with cubic interpolation
__global__
void boundary_location_backup(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
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

// make corrections to xpr etc with 5th order interpolation
// seems that 6th order interpolation gives the same results as cubic interpolation
__global__
void boundary_location(double * xpr, double * xpl, double * ypf, double * ypb, double * zpu, double * zpd, double const * lsf, int num_ele, int rows, int cols, int pges, double dx, double dy, double dz)
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

