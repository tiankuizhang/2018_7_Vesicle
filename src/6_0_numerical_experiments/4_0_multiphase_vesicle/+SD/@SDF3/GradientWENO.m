function [fx,fy,fz] = GradientWENO(obj,F)
	
	% calculate the extend velocity upwindly
	fx = zeros(obj.GD3.Size, 'gpuArray');
	fy = zeros(obj.GD3.Size, 'gpuArray');
	fz = zeros(obj.GD3.Size, 'gpuArray');
	
	[fx, fy, fz] = feval(obj.WENO5RK3_upwind_normal, ...
			fx, fy, fz, obj.F, obj.xpr, obj.xpl, obj.ypf, obj.ypb, obj.zpu, obj.zpd, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);
	
	fgradient = max(sqrt(fx.^2+fy.^2+fz.^2), 1e-6);
	
	nx = fx ./ fgradient;
	ny = fy ./ fgradient;
	nz = fz ./ fgradient;
	
	vx = sign(obj.F) .* nx;
	vy = sign(obj.F) .* ny;
	vz = sign(obj.F) .* nz;

	% interpolate values for crossing points
	cpr = zeros(obj.GD3.Size, 'gpuArray');
	cpl = zeros(obj.GD3.Size, 'gpuArray');
	cpf = zeros(obj.GD3.Size, 'gpuArray');
	cpb = zeros(obj.GD3.Size, 'gpuArray');
	cpu = zeros(obj.GD3.Size, 'gpuArray');
	cpd = zeros(obj.GD3.Size, 'gpuArray');

	[cpr,cpl,cpf,cpb,cpu,cpd] = feval(obj.WENO5RK3_boundary_interpolant, ...
			cpr,cpl,cpf,cpb,cpu,cpd,obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,F, ...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);

	% calculate gradient of field F
	[fx, fy, fz] = feval(obj.upwind_derivative, ...
			fx,fy,fz,F,vx,vy,vz,...
			obj.xpr,obj.xpl,obj.ypf,obj.ypb,obj.zpu,obj.zpd,...
			cpr,cpl,cpf,cpb,cpu,cpd,...
			obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
			obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);

end
