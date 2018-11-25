/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the set calculus toolbox scheme, i.e., to calculate gradients,normal vectors,
 * curvatures, Heaviside function and Dirac_Delta function 
 ******************************************************************************/
#include "shared_utilities.cuh"
#include "shared_utilities.cup"

__global__ 
void set_calculus_toolbox(double * Fx, double * Fy, double * Fz, double * FGradMag,	double * Nx, double * Ny, double * Nz, double * Fxx, double * Fyy, double * Fzz, double * Fxy, double * Fyz, double * Fzx, double * FLaplacian, double * MeanCurvature, double * GaussianCurvature, double * Heaviside, double * DiracDelta, double const * lsf, double const * HPrimal, int rows, int cols, int pges, double dx, double dy, double dz, double ds, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	int right 	= sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int left 	= sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);

	int front 	= sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);	
	int back 	= sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);

	int up 		= sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);	
	int down 	= sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);

	double fx = (lsf[right] - lsf[left]) / (2*dx);
	double fy = (lsf[front] - lsf[back]) / (2*dy);
	double fz = (lsf[up] - lsf[down]) / (2*dz);
	double fGradMag = sqrt(fx*fx + fy*fy + fz*fz); 
	fGradMag = max2(fGradMag, 1e-14); // avoid signularity
	
	Fx[ind] = fx;
	Fy[ind] = fy;
	Fz[ind] = fz;
	FGradMag[ind] = fGradMag;

	Nx[ind] = fx / fGradMag;
	Ny[ind] = fy / fGradMag;
	Nz[ind] = fz / fGradMag;

	int front_right = sub2ind(row_idx+1, col_idx+1, pge_idx, rows, cols, pges);
	int back_left 	= sub2ind(row_idx-1, col_idx-1, pge_idx, rows, cols, pges);
	int back_right 	= sub2ind(row_idx-1, col_idx+1, pge_idx, rows, cols, pges);
	int front_left 	= sub2ind(row_idx+1, col_idx-1, pge_idx, rows, cols, pges);

	int front_up 	= sub2ind(row_idx+1, col_idx, pge_idx+1, rows, cols, pges);
	int back_down 	= sub2ind(row_idx-1, col_idx, pge_idx-1, rows, cols, pges);
	int front_down 	= sub2ind(row_idx+1, col_idx, pge_idx-1, rows, cols, pges);
	int back_up 	= sub2ind(row_idx-1, col_idx, pge_idx+1, rows, cols, pges);

	int right_up 	= sub2ind(row_idx, col_idx+1, pge_idx+1, rows, cols, pges);
	int left_down 	= sub2ind(row_idx, col_idx-1, pge_idx-1, rows, cols, pges);
	int right_down 	= sub2ind(row_idx, col_idx+1, pge_idx-1, rows, cols, pges);
	int left_up 	= sub2ind(row_idx, col_idx-1, pge_idx+1, rows, cols, pges);

	double fxx = (lsf[right] - 2*lsf[ind] + lsf[left]) / (dx*dx);
	double fyy = (lsf[front] - 2*lsf[ind] + lsf[back]) / (dy*dy);
	double fzz = (lsf[up] - 2*lsf[ind] + lsf[down]) / (dz*dz);
	double fxy = (lsf[front_right]+lsf[back_left]-lsf[front_left]-lsf[back_right]) / (4*ds*ds);
	double fyz = (lsf[front_up]+lsf[back_down]-lsf[front_down]-lsf[back_up]) / (4*ds*ds);
	double fzx = (lsf[right_up]+lsf[left_down]-lsf[right_down]-lsf[left_up]) / (4*ds*ds);
	double fLaplacian = fxx + fyy + fzz;

	Fxx[ind] = fxx;
	Fyy[ind] = fyy;
	Fzz[ind] = fzz;
	Fxy[ind] = fxy;
	Fyz[ind] = fyz;
	Fzx[ind] = fzx;
	FLaplacian[ind] = fLaplacian;

	// calculate mean curvature
	double col1 = fxx*fx + fxy*fy + fzx*fz;
	double col2 = fxy*fx + fyy*fy + fyz*fz;
	double col3 = fzx*fx + fyz*fy + fzz*fz;
	
	MeanCurvature[ind] = - fLaplacian/fGradMag + (fx*col1+fy*col2+fz*col3)/pow(fGradMag,3);

	// calculate Gaussian curvature
	col1 = (fyy*fzz-fyz*fyz)*fx + (fzx*fyz-fxy*fzz)*fy + (fxy*fyz-fzx*fyy)*fz;
	col2 = (fyz*fzx-fxy*fzz)*fx + (fxx*fzz-fzx*fzx)*fy + (fzx*fxy-fxx*fyz)*fz;
	col3 = (fxy*fyz-fyy*fzx)*fx + (fzx*fxy-fxx*fyz)*fy + (fxx*fyy-fxy*fxy)*fz;

	GaussianCurvature[ind] = (fx*col1+fy*col2+fz*col3) / pow(fGradMag,4);

	// calculate Heaviside function
	double px = (HPrimal[right] - HPrimal[left]) / (2*dx);
	double py = (HPrimal[front] - HPrimal[back]) / (2*dy);
	double pz = (HPrimal[up] - HPrimal[down]) / (2*dz);

	double dot_DHPrimal_DF = px*fx + py*fy + pz*fz;

	Heaviside[ind] = dot_DHPrimal_DF / pow(fGradMag,2);

	// calculate DiraDelta function
	double pxx = (HPrimal[right] - 2*HPrimal[ind] +HPrimal[left]) / (dx*dx);
	double pyy = (HPrimal[front] - 2*HPrimal[ind] + HPrimal[back]) / (dy*dy);
	double pzz = (HPrimal[up] - 2*HPrimal[ind] + HPrimal[down]) / (dz*dz);
	double pLaplacian = pxx + pyy + pzz;

	DiracDelta[ind] = pLaplacian/pow(fGradMag,2) - dot_DHPrimal_DF*fLaplacian/pow(fGradMag,4);
}	



















