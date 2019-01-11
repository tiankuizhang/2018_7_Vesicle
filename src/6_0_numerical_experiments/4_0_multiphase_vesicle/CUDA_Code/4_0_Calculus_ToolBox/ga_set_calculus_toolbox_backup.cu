/*********************************************************************************
 * ga_set_calculus_toolbox.cu // ga stands for geometry aware
 * calculate gradient, hessian, curvature of the level set function
 * with 4th order central scheme if there is no kink within the stencil
 * otherwise, use 2nd order central scheme
 * note that it seems very important to use a compact stencil 
 ********************************************************************************/

#include "shared_utilities.cuh"
#include "shared_utilities.cup"

__global__
void ga_set_calculus_toolbox(double * Fx, double * Fy, double * Fz, double * FGradMag, double * Nx, double * Ny, double * Nz, double * Fxx, double * Fyy, double * Fzz, double * Fxy, double * Fyz, double * Fzx, double * FLaplacian, double * MeanCurvature, double * GaussianCurvature, double * Heaviside, double * DiracDelta, double const * lsf, bool const * kink, double const * HPrimal, int rows, int cols, int pges, double dx, double dy, double dz, double ds, int num_ele)
{
	int row_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int col_idx = blockIdx.y * blockDim.y + threadIdx.y;
	int pge_idx = blockIdx.z * blockDim.z + threadIdx.z;

	if(row_idx >= rows || col_idx >= cols || pge_idx >= pges){
		return;
	}

	int ind = sub2ind(row_idx, col_idx, pge_idx, rows, cols, pges);

	bool local_kink;

	int rght1 = sub2ind(row_idx, col_idx+1, pge_idx, rows, cols, pges);
	int rght2 = sub2ind(row_idx, col_idx+2, pge_idx, rows, cols, pges);
	int left1 = sub2ind(row_idx, col_idx-1, pge_idx, rows, cols, pges);
	int left2 = sub2ind(row_idx, col_idx-2, pge_idx, rows, cols, pges);

	int frnt1 = sub2ind(row_idx+1, col_idx, pge_idx, rows, cols, pges);
	int frnt2 = sub2ind(row_idx+2, col_idx, pge_idx, rows, cols, pges);
	int back1 = sub2ind(row_idx-1, col_idx, pge_idx, rows, cols, pges);
	int back2 = sub2ind(row_idx-2, col_idx, pge_idx, rows, cols, pges);

	int upup1 = sub2ind(row_idx, col_idx, pge_idx+1, rows, cols, pges);
	int upup2 = sub2ind(row_idx, col_idx, pge_idx+2, rows, cols, pges);
	int down1 = sub2ind(row_idx, col_idx, pge_idx-1, rows, cols, pges);
	int down2 = sub2ind(row_idx, col_idx, pge_idx-2, rows, cols, pges);

	double fx2 = (lsf[rght1] - lsf[left1]) / (2.0*dx);
	double fy2 = (lsf[frnt1] - lsf[back1]) / (2.0*dy);
	double fz2 = (lsf[upup1] - lsf[down1]) / (2.0*dz);

	double fx4 = (-lsf[rght2] + 8.0*lsf[rght1] - 8.0*lsf[left1] + lsf[left2]) / (12.0*dx);
	double fy4 = (-lsf[frnt2] + 8.0*lsf[frnt1] - 8.0*lsf[back1] + lsf[back2]) / (12.0*dy);
	double fz4 = (-lsf[upup2] + 8.0*lsf[upup1] - 8.0*lsf[down1] + lsf[down2]) / (12.0*dz);

	double fxx2 = (lsf[rght1] - 2.0*lsf[ind] + lsf[left1]) / (dx*dx);
	double fyy2 = (lsf[frnt1] - 2.0*lsf[ind] + lsf[back1]) / (dy*dy);
	double fzz2 = (lsf[upup1] - 2.0*lsf[ind] + lsf[down1]) / (dz*dz);

	double fxx4 = (-lsf[rght2] + 16.0*lsf[rght1] - 30.0 * lsf[ind] + 16.0*lsf[left1] - lsf[left2]) / (12.0*dx*dx);
	double fyy4 = (-lsf[frnt2] + 16.0*lsf[frnt1] - 30.0 * lsf[ind] + 16.0*lsf[back1] - lsf[back2]) / (12.0*dy*dy);
	double fzz4 = (-lsf[upup2] + 16.0*lsf[upup1] - 30.0 * lsf[ind] + 16.0*lsf[down1] - lsf[down2]) / (12.0*dz*dz);

	local_kink = kink[rght2] || kink[rght1] || kink[ind] || kink[left1] || kink[left2]; 
	double fx = local_kink ? fx2 : fx4;
	double fxx = local_kink ? fxx2 : fxx4;

	local_kink = kink[frnt2] || kink[frnt1] || kink[ind] || kink[back1] || kink[back2];
	double fy = local_kink ? fy2 : fy4;
	double fyy = local_kink ? fyy2 : fyy4;

	local_kink = kink[upup2] || kink[upup1] || kink[ind] || kink[down1] || kink[down2];
	double fz = local_kink ? fz2 : fz4;
	double fzz = local_kink ? fzz2 : fzz4;

	double fLaplacian = fxx + fyy + fzz;

	double fGradMag = sqrt(fx*fx + fy*fy + fz*fz); 
	fGradMag = max2(fGradMag, 1e-14); // avoid signularity

	Fx[ind] = fx;
	Fy[ind] = fy;
	Fz[ind] = fz;
	Fxx[ind] = fxx;
	Fyy[ind] = fyy;
	Fzz[ind] = fzz;
	FLaplacian[ind] = fLaplacian;
	FGradMag[ind] = fGradMag;

	Nx[ind] = fx / fGradMag;
	Ny[ind] = fy / fGradMag;
	Nz[ind] = fz / fGradMag;

	int frnt_rght1 = sub2ind(row_idx+1, col_idx+1, pge_idx, rows, cols, pges);
	int frnt_rght2 = sub2ind(row_idx+2, col_idx+2, pge_idx, rows, cols, pges);
	int back_left1 = sub2ind(row_idx-1, col_idx-1, pge_idx, rows, cols, pges);
	int back_left2 = sub2ind(row_idx-2, col_idx-2, pge_idx, rows, cols, pges);
	int back_rght1 = sub2ind(row_idx-1, col_idx+1, pge_idx, rows, cols, pges);
	int back_rght2 = sub2ind(row_idx-2, col_idx+2, pge_idx, rows, cols, pges);
	int frnt_left1 = sub2ind(row_idx+1, col_idx-1, pge_idx, rows, cols, pges);
	int frnt_left2 = sub2ind(row_idx+2, col_idx-2, pge_idx, rows, cols, pges);

	int frnt_upup1 = sub2ind(row_idx+1, col_idx, pge_idx+1, rows, cols, pges);
	int frnt_upup2 = sub2ind(row_idx+2, col_idx, pge_idx+2, rows, cols, pges);
	int back_down1 = sub2ind(row_idx-1, col_idx, pge_idx-1, rows, cols, pges);
	int back_down2 = sub2ind(row_idx-2, col_idx, pge_idx-2, rows, cols, pges);
	int frnt_down1 = sub2ind(row_idx+1, col_idx, pge_idx-1, rows, cols, pges);
	int frnt_down2 = sub2ind(row_idx+2, col_idx, pge_idx-2, rows, cols, pges);
	int back_upup1 = sub2ind(row_idx-1, col_idx, pge_idx+1, rows, cols, pges);
	int back_upup2 = sub2ind(row_idx-2, col_idx, pge_idx+2, rows, cols, pges);

	int rght_upup1 = sub2ind(row_idx, col_idx+1, pge_idx+1, rows, cols, pges);
	int rght_upup2 = sub2ind(row_idx, col_idx+2, pge_idx+2, rows, cols, pges);
	int left_down1 = sub2ind(row_idx, col_idx-1, pge_idx-1, rows, cols, pges);
	int left_down2 = sub2ind(row_idx, col_idx-2, pge_idx-2, rows, cols, pges);
	int rght_down1 = sub2ind(row_idx, col_idx+1, pge_idx-1, rows, cols, pges);
	int rght_down2 = sub2ind(row_idx, col_idx+2, pge_idx-2, rows, cols, pges);
	int left_upup1 = sub2ind(row_idx, col_idx-1, pge_idx+1, rows, cols, pges);
	int left_upup2 = sub2ind(row_idx, col_idx-2, pge_idx+2, rows, cols, pges);

	double fxy2 = (lsf[frnt_rght1]+lsf[back_left1]-lsf[frnt_left1]-lsf[back_rght1]) / (4*ds*ds);
	double fyz2 = (lsf[frnt_upup1]+lsf[back_down1]-lsf[frnt_down1]-lsf[back_upup1]) / (4*ds*ds);
	double fzx2 = (lsf[rght_upup1]+lsf[left_down1]-lsf[rght_down1]-lsf[left_upup1]) / (4*ds*ds);

	double fxy4 = (-lsf[frnt_rght2]-lsf[back_left2]+lsf[frnt_left2]+lsf[back_rght2]+16.0*lsf[frnt_rght1]+16.0*lsf[back_left1]-16.0*lsf[frnt_left1]-16.0*lsf[back_rght1]) / (48.0*dx*dy);
	double fyz4 = (-lsf[frnt_upup2]-lsf[back_down2]+lsf[frnt_down2]+lsf[back_upup2]+16.0*lsf[frnt_upup1]+16.0*lsf[back_down1]-16.0*lsf[frnt_down1]-16.0*lsf[back_upup1]) / (48.0*dy*dz);
	double fzx4 = (-lsf[rght_upup2]-lsf[left_down2]+lsf[rght_down2]+lsf[left_upup2]+16.0*lsf[rght_upup1]+16.0*lsf[left_down1]-16.0*lsf[rght_down1]-16.0*lsf[left_upup1]) / (48.0*dz*dx);

	local_kink = kink[frnt_rght2] || kink[back_left2] || kink[frnt_left2] || kink[back_rght2] || kink[frnt_rght1] || kink[back_left1] || kink[frnt_left1] || kink[back_rght1] || kink[ind];
	double fxy = local_kink ? fxy2 : fxy4;

	local_kink = kink[frnt_upup2] || kink[back_down1] || kink[frnt_down2] || kink[back_upup2] || kink[frnt_upup1] || kink[back_down1] || kink[frnt_down1] || kink[back_upup1] || kink[ind];
	double fyz = local_kink ? fyz2 : fyz4;

	local_kink = kink[rght_upup2] || kink[left_down2] || kink[rght_down2] || kink[left_upup2] || kink[rght_upup1] || kink[left_down1] || kink[rght_down1] || kink[left_upup1] || kink[ind];
	double fzx = local_kink ? fzx2 : fzx4;

	Fxy[ind] = fxy;
	Fyz[ind] = fyz;
	Fzx[ind] = fzx;

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
	double px = (HPrimal[rght1] - HPrimal[left1]) / (2*dx);
	double py = (HPrimal[frnt1] - HPrimal[back1]) / (2*dy);
	double pz = (HPrimal[upup1] - HPrimal[down1]) / (2*dz);

	double dot_DHPrimal_DF = px*fx + py*fy + pz*fz;

	Heaviside[ind] = dot_DHPrimal_DF / pow(fGradMag,2);

	// calculate DiraDelta function
	double pxx = (HPrimal[rght1] - 2*HPrimal[ind] + HPrimal[left1]) / (dx*dx);
	double pyy = (HPrimal[frnt1] - 2*HPrimal[ind] + HPrimal[back1]) / (dy*dy);
	double pzz = (HPrimal[upup1] - 2*HPrimal[ind] + HPrimal[down1]) / (dz*dz);
	double pLaplacian = pxx + pyy + pzz;

	DiracDelta[ind] = pLaplacian/pow(fGradMag,2) - dot_DHPrimal_DF*fLaplacian/pow(fGradMag,4);

}
