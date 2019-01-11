/*******************************************************************************
 * serveral useful gpu functions will be defined in this file to facilitate
 * the set calculus toolbox scheme, i.e., to calculate gradients,normal vectors,
 * curvatures, Heaviside function and Dirac_Delta function 
 ******************************************************************************/
#include "shared_utilities.cuh"
#include "shared_utilities.cup"


__global__
void auxi_set_calculus_toolbox(double * Ax, double * Ay, double * Az, double * AGradMag, double * ACx, double * ACy, double * ACz, double * ANormCrossAF, double * Tx, double * Ty, double * Tz, double * Anx, double * Any, double * Anz, double * Axx, double * Ayy, double * Azz, double * Axy, double * Ayz, double * Azx, double * ALaplacian, double * GeodesicCurvature, double * NormalCurvature, double * GeodesicTorsion, double * BPerpendicular, double * AHeaviside, double * ADiracDelta, double const * lsf, double const * AHPrimal, double const * Fx, double const * Fy, double const * Fz, double const * FGradMag, double const * Nx, double const * Ny, double const * Nz, double const * Fxx, double const * Fyy, double const * Fzz, double const * Fxy, double const * Fyz, double const * Fzx, int rows, int cols, int pges, double dx, double dy, double dz, double ds, int num_ele)
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

	double ax = (lsf[right] - lsf[left]) / (2*dx);
	double ay = (lsf[front] - lsf[back]) / (2*dy);
	double az = (lsf[up] - lsf[down]) / (2*dz);
	double aGradMag = norm(ax, ay, az);

	Ax[ind] = ax;
	Ay[ind] = ay;
	Az[ind] = az;
	AGradMag[ind] = aGradMag;

	double fx = Fx[ind];
	double fy = Fy[ind];
	double fz = Fz[ind];

	double Cx, Cy, Cz;
	cross_product(Cx,Cy,Cz,fx,fy,fz,ax,ay,az);
	double NormCrossAF = norm(Cx,Cy,Cz);

	ACx[ind] = Cx;
	ACy[ind] = Cy;
	ACz[ind] = Cz;
	ANormCrossAF[ind] = NormCrossAF;

	double tx = Cx / NormCrossAF;
	double ty = Cy / NormCrossAF;
	double tz = Cz / NormCrossAF;

	Tx[ind] = tx;
	Ty[ind] = ty;
	Tz[ind] = tz;

	double fNx = Nx[ind];
	double fNy = Ny[ind];
	double fNz = Nz[ind];

	double nx, ny, nz;
	cross_product(nx,ny,nz,tx,ty,tz,fNx,fNy,fNz);

	Anx[ind] = nx;
	Any[ind] = ny;
	Anz[ind] = nz;

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

	double axx = (lsf[right] - 2*lsf[ind] + lsf[left]) / (dx*dx);
	double ayy = (lsf[front] - 2*lsf[ind] + lsf[back]) / (dy*dy);
	double azz = (lsf[up] - 2*lsf[ind] + lsf[down]) / (dz*dz);
	double axy = (lsf[front_right]+lsf[back_left]-lsf[front_left]-lsf[back_right]) / (4*ds*ds);
	double ayz = (lsf[front_up]+lsf[back_down]-lsf[front_down]-lsf[back_up]) / (4*ds*ds);
	double azx = (lsf[right_up]+lsf[left_down]-lsf[right_down]-lsf[left_up]) / (4*ds*ds);
	double aLaplacian = axx + ayy + azz;

	Axx[ind] = axx;
	Ayy[ind] = ayy;
	Azz[ind] = azz;
	Axy[ind] = axy;
	Ayz[ind] = ayz;
	Azx[ind] = azx;
	ALaplacian[ind] = aLaplacian;

	// geodesic curvature
	double fxx = Fxx[ind];
	double fyy = Fyy[ind];
	double fzz = Fzz[ind];
	double fxy = Fxy[ind];
	double fyz = Fyz[ind];
	double fzx = Fzx[ind];
	double fGradMag = FGradMag[ind];

	double vx = tx*fxx + ty*fxy + tz*fzx;
	double vy = tx*fxy + ty*fyy + tz*fyz;
	double vz = tx*fzx + ty*fyz + tz*fzz;

	double w1x, w1y, w1z;
	cross_product(w1x,w1y,w1z,vx,vy,vz,ax,ay,az);

	vx = tx*axx + ty*axy + tz*azx;
	vy = tx*axy + ty*ayy + tz*ayz;
	vz = tx*azx + ty*ayz + tz*azz;

	double w2x, w2y, w2z;
	cross_product(w2x,w2y,w2z,fx,fy,fz,vx,vy,vz);

	GeodesicCurvature[ind] = ( nx*(w1x+w2x) + ny*(w1y+w2y) + nz*(w1z+w2z) ) / NormCrossAF;

	/* NormalCurvature, GeodesicTorsion, BPerpendicular */
	double Nxx = fxx / fGradMag - fx*(fxx*fx + fxy*fy + fzx*fz) / pow(fGradMag,3) ; 
	double Nyx = fxy / fGradMag - fy*(fxx*fx + fxy*fy + fzx*fz) / pow(fGradMag,3) ; 
	double Nzx = fzx / fGradMag - fz*(fxx*fx + fxy*fy + fzx*fz) / pow(fGradMag,3) ; 

	double Nxy = fxy / fGradMag - fx*(fxy*fx + fyy*fy + fyz*fz) / pow(fGradMag,3) ;
	double Nyy = fyy / fGradMag - fy*(fxy*fx + fyy*fy + fyz*fz) / pow(fGradMag,3) ;
	double Nzy = fyz / fGradMag - fz*(fxy*fx + fyy*fy + fyz*fz) / pow(fGradMag,3) ;

	double Nxz = fzx / fGradMag - fx*(fzx*fx + fyz*fy + fzz*fz) / pow(fGradMag,3) ;
	double Nyz = fyz / fGradMag - fy*(fzx*fx + fyz*fy + fzz*fz) / pow(fGradMag,3) ;
	double Nzz = fzz / fGradMag - fz*(fzx*fx + fyz*fy + fzz*fz) / pow(fGradMag,3) ;

	// NormalCurvature.
	vx = Nxx * tx + Nxy * ty + Nxz * tz; 
	vy = Nyx * tx + Nyy * ty + Nyz * tz; 
	vz = Nzx * tx + Nzy * ty + Nzz * tz; 

	NormalCurvature[ind] = - (tx*vx + ty*vy + tz*vz);

	// GeodesicTorsion, BPerpendicular
	vx = Nxx * nx + Nxy * ny + Nxz * nz; 
	vy = Nyx * nx + Nyy * ny + Nyz * nz; 
	vz = Nzx * nx + Nzy * ny + Nzz * nz; 

	GeodesicTorsion[ind] = - (tx*vx + ty*vy + tz*vz);
	BPerpendicular[ind] = - (nx*vx + ny*vy + nz*vz);

	/*primal of Heaviside(A), Heaviside(A), DiracDelta(A)*/

	// calculate Heaviside function
	double px = (AHPrimal[right] - AHPrimal[left]) / (2*dx);
	double py = (AHPrimal[front] - AHPrimal[back]) / (2*dy);
	double pz = (AHPrimal[up] - AHPrimal[down]) / (2*dz);

	double dot_DAHPrimal_DF = px*ax + py*ay + pz*az;

	AHeaviside[ind] = dot_DAHPrimal_DF / pow(aGradMag,2);

	// calculate DiraDelta function
	double pxx = (AHPrimal[right] - 2*AHPrimal[ind] +AHPrimal[left]) / (dx*dx);
	double pyy = (AHPrimal[front] - 2*AHPrimal[ind] + AHPrimal[back]) / (dy*dy);
	double pzz = (AHPrimal[up] - 2*AHPrimal[ind] + AHPrimal[down]) / (dz*dz);
	double pLaplacian = pxx + pyy + pzz;

	ADiracDelta[ind] = pLaplacian/pow(aGradMag,2) - dot_DAHPrimal_DF*aLaplacian/pow(aGradMag,4);

}






















