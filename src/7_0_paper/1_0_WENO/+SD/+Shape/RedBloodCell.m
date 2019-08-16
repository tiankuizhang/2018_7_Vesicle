% return shape of a red blood cell
function [X,Y,Z,f,Nx,Ny,Nz,MC,GC,SL_MC] = RedBloodCell(Size,R)
	
	% create a meshgrid
	Nx = Size(1);
	Ny = Size(2);
	Nz = Size(3);

	xmin = -1.3*R;
	xmax =  1.3*R;
	xv = linspace(xmin,xmax,Nx);
	dx = xv(2) - xv(1);

	yv = ( (1-Ny)/2. : (Ny-1)/2 ) * dx;
	zv = ( (1-Nz)/2. : (Nz-1)/2 ) * dx;

	[X,Y,Z] = meshgrid(xv,yv,zv);

	X = gpuArray(X);
	Y = gpuArray(Y);
	Z = gpuArray(Z);

	% create the level set function representing the red blood cell
	syms x y z
	C0 = 0.2072;
	C1 = 2.0026;
	C2 = -1.1228;

	f = 4*z.^2/R^2 - (1-(x.^2+y.^2)/R^2) .* (C0 + C1*(x.^2+y.^2)/R^2 + C2*(x.^2+y.^2).^2/R^4).^2;

	H = hessian(f,[x,y,z]);
	G = gradient(f,[x,y,z]);

	fxx = H(1,1); fxy = H(1,2); fxz = H(1,3);
	fyx = H(2,1); fyy = H(2,2); fyz = H(2,3);
	fzx = H(3,1); fzy = H(3,2); fzz = H(3,3);

	% cofactor matrix of Hessian
	cH = [fyy*fzz-fzy*fyz, fzx*fyz-fyx*fzz, fyx*fzy-fzx*fyy;
		  fzy*fxz-fxy*fzz, fxx*fzz-fzx*fxz, fzx*fxy-fxx*fzy;
		  fxy*fyz-fyy*fxz, fxz*fyx-fxx*fyz, fxx*fyy-fyx*fxy ];

	% normals
	FGradMag = sqrt(G'*G+eps);
	N = G / FGradMag;

	% mean curvature
	MC = -trace(H) / FGradMag + G' * H * G / FGradMag^3;

	% surface laplacian of mean curvature
	G_MC = gradient(MC,[x,y,z]);
	H_MC = hessian(MC,[x,y,z]);
	SL_MC = trace(H_MC) - N' * H_MC * N - N' * G_MC * MC; 

	% Gaussian curvature
	GC = G' * cH * G / FGradMag^4;

	f = matlabFunction(f);
	Nx = matlabFunction(N(1));
	Ny = matlabFunction(N(2));
	Nz = matlabFunction(N(3));
	MC = matlabFunction(MC);
	GC = matlabFunction(GC);
	SL_MC = matlabFunction(SL_MC);

end

