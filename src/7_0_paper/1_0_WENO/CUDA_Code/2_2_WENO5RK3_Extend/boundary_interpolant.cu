/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the extension scheme 
 ******************************************************************************/
#include "shared_utilities.cuh"
#include "shared_utilities.cup"

// calculate the WENO interpolant at cx, given data (x_m2,y_m2),...,(x3,y3)
// aasuem x0 < cx < x1
__device__ inline
void weno_interpolant(double & c, double cx, double s, double const * x, double const * y)
{
	double epsilon = 0.000001;
	double xm2 = x[0], xm1 = x[1], x0 = x[2], x1 = x[3], x2 = x[4], x3 = x[5]; 
	double ym2 = y[0], ym1 = y[1], y0 = y[2], y1 = y[3], y2 = y[4], y3 = y[5]; 

	double p1, p2, p3, c1, c2, c3, IS1, IS2, IS3; // ENO candidates

	p1 = ym2 	+ (ym1 - ym2) * (cx - xm2) / s 
				+ (y0 - 2*ym1 + ym2) * (cx - xm2) * (cx - xm1) / (2.0 * s * s)
				+ (y1 - 3*y0 + 3*ym1 - ym2) * (cx - xm2) * (cx - xm1) * (cx - x0) / (6 * s * s * s);

	p2 = ym1 	+ (y0 - ym1) * (cx - xm1) / s 
				+ (y1 - 2*y0 + ym1) * (cx - xm1) * (cx - x0) / (2.0 * s * s)
				+ (y2 - 3*y1 + 3*y0 - ym1) * (cx - xm1) * (cx - x0) * (cx - x1) / (6 * s * s * s);

	p3 = y0 	+ (y1 - y0) * (cx - x0) / s 
				+ (y2 - 2*y1 + y0) * (cx - x0) * (cx - x1) / (2.0 * s * s)
				+ (y3 - 3*y2 + 3*y1 - y0) * (cx - x1) * (cx - x1) * (cx - x2) / (6 * s * s * s);

	c1 = (x2 - cx) * (x3 - cx)   / (20 * s * s);
	c2 = (x3 - cx) * (cx - xm2)  / (10 * s * s);
	c3 = (cx - xm2) * (cx - xm1) / (20 * s * s);

	IS1 = (-3579*y1*y0 + 2634*y1*ym1 - 683*y1*ym2 - 6927*y0*ym1 + 1854*y0*ym2 - 1659*ym1*ym2 + 814*y1*y1 + 4326*y0*y0 + 2976*ym1*ym1 + 244*ym2*ym2) / 180;
	IS2 = (-3777*y1*y0 + 1074*y1*ym1 + 1269*y0*ym1 + 1986*y1*y1 + 1986*y0*y0 + 244*ym1*ym1 + 244*y2*y2 - 1269*y2*y1 + 1074*y2*y0 - 293*y2*ym1) / 180;
	IS3 = (-3579*y1*y0 + 4326*y1*y1 + 814*y0*y0 + 2976*y2*y2 + 244*y3*y3 - 683*y3*y0 - 6927*y2*y1 + 2634*y2*y0 - 1659*y3*y2 + 1854*y3*y1) / 180;

	double alpha1, alpha2, alpha3, sum, omega1, omega2, omega3;
	alpha1 = c1 / pow( (IS1 + epsilon), 2);
	alpha2 = c2 / pow( (IS2 + epsilon), 2);
	alpha3 = c3 / pow( (IS3 + epsilon), 2);
	sum = alpha1 + alpha2 + alpha3;
	omega1 = alpha1 / sum;
	omega2 = alpha2 / sum;
	omega3 = alpha3 / sum;

	c = omega1 * p1 + omega2 * p2 + omega3 * p3;

}

// suppose that (3*s,v3) is a boundary node and we know the distance from this node to the boundary
// sixth_interp will interpolate values onto the boundary node from the stencil
// (0,v0),(s,v1),(2s,v2),(3s,v3),(4s,v4),(5s,v5) if p2*p3<0 for c_backward
// (s,v1),(2s,v2),(3s,v3),(4s,v4),(5s,v5),(6s,v6) if p3*p4<0 for c_forward
__device__ inline
void sixth_interp_backup(double & c_forward, double & c_backward, double dis_f, double dis_b, double s, double v0, double v1, double v2, double v3, double v4, double v5, double v6)
{
	c_forward = 0.0;
	c_backward = 0.0;

	// if there is a boundary in the forward direction
	if(dis_f<s){
		double xf[6] = {s,2*s,3*s,4*s,5*s,6*s};
		double yf[6] = {v1,v2,v3,v4,v5,v6};
		weno_interpolant(c_forward,3*s+dis_f,s,xf,yf);
	}
	// if there is a boundary in the backward direction
	if(dis_b<s){
		double xf[6] = {0.0,s,2*s,3*s,4*s,5*s};
		double yf[6] = {v0,v1,v2,v3,v4,v5};
		double cx = 3*s - dis_b;
		weno_interpolant(c_backward,3*s-dis_b,s,xf,yf);
	}
}

__device__ inline
void sixth_interp(double & c_forward, double & c_backward, double dis_f, double dis_b, double s, double v0, double v1, double v2, double v3, double v4, double v5, double v6)
{
	c_forward = 0.0;
	c_backward = 0.0;

	// if there is a boundary in the forward direction
	if(dis_f<s){
		double xf[6] = {s,2*s,3*s,4*s,5*s,6*s};
		double yf[6] = {v1,v2,v3,v4,v5,v6};
		Poly<5> p5_f;
		p5_f.setInterpolation(xf,yf);
		c_forward = p5_f(3*s + dis_f);
	}
	// if there is a boundary in the backward direction
	if(dis_b<s){
		double xf[6] = {0.0,s,2*s,3*s,4*s,5*s};
		double yf[6] = {v0,v1,v2,v3,v4,v5};
		Poly<5> p5_b;
		p5_b.setInterpolation(xf,yf);
		c_backward = p5_b(3*s - dis_b);
	}
}

// interpolate values at boundary points
__global__
void boundary_interpolant(double * cpr, double * cpl, double * cpf, double * cpb, double * cpu, double * cpd, double const * xpr, double const * xpl, double const * ypf, double const * ypb, double const * zpu, double const * zpd, double const * lsf, int rows, int cols, int pges, double dx, double dy, double dz, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int right  = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int right2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int right3 = sub2ind(row_idx, col_idx+3, pge_idx, rows, cols, pges);
	int left   = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2  = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);
	int left3  = sub2ind(row_idx, col_idx-3, pge_idx, rows, cols, pges);

	sixth_interp(cpr[ind],cpl[ind],xpr[ind],xpl[ind],dx,lsf[left3],lsf[left2],lsf[left],lsf[ind],lsf[right],lsf[right2],lsf[right3]);

	int front  = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	int front2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int front3 = sub2ind(row_idx+3, col_idx, pge_idx, rows, cols, pges);
	int back   = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2  = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);
	int back3  = sub2ind(row_idx-3, col_idx, pge_idx, rows, cols, pges);

	sixth_interp(cpf[ind],cpb[ind],ypf[ind],ypb[ind],dy,lsf[back3],lsf[back2],lsf[back],lsf[ind],lsf[front],lsf[front2],lsf[front3]);

	int up    = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	int up2   = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
	int up3   = sub2ind(row_idx, col_idx, pge_idx+3, rows, cols, pges);
	int down  = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);
	int down3 = sub2ind(row_idx, col_idx, pge_idx-3, rows, cols, pges);

	sixth_interp(cpu[ind],cpd[ind],zpu[ind],zpd[ind],dz,lsf[down3],lsf[down2],lsf[down],lsf[ind],lsf[up],lsf[up2],lsf[up3]);
}







