% test FFT scheme

n = 128;

xv = linspace(-2,2,n);
yv = xv;
zv = xv;

[x,y,z] = meshgrid(xv,yv,zv);

x = gpuArray(x);
y = gpuArray(y);
z = gpuArray(z);

grid = SD.GD3(x,y,z);

F = x.^2 + 2 * y.^2;

a = - 0.01;
b = 0.01;

%Opgpu = a * grid.LLaplacian + grid.Idt;
Opgpu = b * grid.LBiLaplacian + grid.Idt;
Opcpu = gather(Opgpu);

%F = rand(n^3,1,'gpuArray');
FCol = reshape(F, [n^3, 1]);

S = Opgpu * FCol;
Scpu = gather(S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test iterative methods

%tic
%[L, U] = ilu(Opcpu, struct('type','nofill','milu','row'));
%toc
%
%tic
%F1 = gmres(@(x)gather(Opgpu*x), Scpu, 10, 1e-6, 300, L, U);
%toc
%F1Col = reshape(F1, grid.Size);
%
%tic
%F2 = bicgstab(@(x)gather(Opgpu*x), Scpu, 1e-6, 200, L, U);
%toc
%F2Col = reshape(F2, grid.Size);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use FFT to solve it

tic
fftS = fftn(reshape(S,grid.Size));
%fftS = fftS ./ (1 - a * ( kx.^2 + ky.^2 + kz.^2) );
fftS = fftS ./ (1 + b * ( grid.kx.^2 + grid.ky.^2 + grid.kz.^2).^2 );
F3 = ifftn(fftS);
F3 = real(F3);
toc

tic
F4 = run_gmres(Opgpu, S, grid, b);
toc

function x1 = run_gmres(Opgpu, Sgpu, grid, b)

	%x1 = gmres(@afun, Sgpu, [], 1e-6, 300, @mfun);
	x1 = gmres(Opgpu, Sgpu, [], 1e-12, 300, @mfun);

	function y = afun(x)
		y = Opgpu * x;
	end

	function y = mfun(S)
		fftS = fftn(reshape(S,grid.Size));
		fftS = fftS ./ (1 + b * ( grid.kx.^2 + grid.ky.^2 + grid.kz.^2).^2 );
		y = real(ifftn(fftS));
		y = reshape(y, [grid.NumElt, 1]);
	end
	
end















