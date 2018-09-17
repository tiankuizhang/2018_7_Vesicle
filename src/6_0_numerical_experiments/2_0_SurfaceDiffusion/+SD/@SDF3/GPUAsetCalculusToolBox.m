% gpu kernel version of AsetCalculusToolBox
function GPUAsetCalculusToolBox(obj)

	obj.AHPrimal = max(obj.A,0);

	obj.Ax = zeros(obj.GD3.Size,'gpuArray');
	obj.Ay = zeros(obj.GD3.Size,'gpuArray');
	obj.Az = zeros(obj.GD3.Size,'gpuArray');
	obj.AGradMag = zeros(obj.GD3.Size,'gpuArray');

	obj.Cx = zeros(obj.GD3.Size,'gpuArray');
	obj.Cy = zeros(obj.GD3.Size,'gpuArray');
	obj.Cz = zeros(obj.GD3.Size,'gpuArray');
	obj.NormCrossAF = zeros(obj.GD3.Size,'gpuArray');

	obj.tx = zeros(obj.GD3.Size,'gpuArray');
	obj.ty = zeros(obj.GD3.Size,'gpuArray');
	obj.tz = zeros(obj.GD3.Size,'gpuArray');

	obj.nx = zeros(obj.GD3.Size,'gpuArray');
	obj.ny = zeros(obj.GD3.Size,'gpuArray');
	obj.nz = zeros(obj.GD3.Size,'gpuArray');

	obj.Axx = zeros(obj.GD3.Size,'gpuArray');
	obj.Ayy = zeros(obj.GD3.Size,'gpuArray');
	obj.Azz = zeros(obj.GD3.Size,'gpuArray');
	obj.Axy = zeros(obj.GD3.Size,'gpuArray');
	obj.Ayz = zeros(obj.GD3.Size,'gpuArray');
	obj.Azx = zeros(obj.GD3.Size,'gpuArray');
	obj.ALaplacian = zeros(obj.GD3.Size,'gpuArray');

	obj.GeodesicCurvature = zeros(obj.GD3.Size,'gpuArray');
	obj.NormalCurvature = zeros(obj.GD3.Size,'gpuArray');
	obj.GeodesicTorsion = zeros(obj.GD3.Size,'gpuArray');
	obj.BPerpendicular = zeros(obj.GD3.Size,'gpuArray');

	obj.AHeaviside = zeros(obj.GD3.Size,'gpuArray');
	obj.ADiracDelta = zeros(obj.GD3.Size,'gpuArray');

	[obj.Ax,obj.Ay,obj.Az,obj.AGradMag,obj.Cx,obj.Cy,obj.Cz,obj.NormCrossAF,...
	 obj.tx,obj.ty,obj.tz,obj.nx,obj.ny,obj.nz,...
	 obj.Axx,obj.Ayy,obj.Azz,obj.Axy,obj.Ayz,obj.Azx,obj.ALaplacian,...
	 obj.GeodesicCurvature,obj.NormalCurvature,obj.GeodesicTorsion,obj.BPerpendicular,...
	 obj.AHeaviside,obj.ADiracDelta] = feval(obj.auxi_set_calculus_toolbox, ...
			obj.Ax,obj.Ay,obj.Az,obj.AGradMag,obj.Cx,obj.Cy,obj.Cz,obj.NormCrossAF,...
	 		obj.tx,obj.ty,obj.tz,obj.nx,obj.ny,obj.nz,...
	 		obj.Axx,obj.Ayy,obj.Azz,obj.Axy,obj.Ayz,obj.Azx,obj.ALaplacian,...
	 		obj.GeodesicCurvature,obj.NormalCurvature,obj.GeodesicTorsion,obj.BPerpendicular,...
			obj.AHeaviside,obj.ADiracDelta,...
			obj.A,obj.AHPrimal,...
			obj.Fx,obj.Fy,obj.Fz,obj.FGradMag,obj.Nx,obj.Ny,obj.Nz,...
			obj.Fxx,obj.Fyy,obj.Fzz,obj.Fxy,obj.Fyz,obj.Fzx,...
			obj.GD3.mrows,obj.GD3.ncols,obj.GD3.lshts,...
			obj.GD3.Dx,obj.GD3.Dy,obj.GD3.Dz,obj.GD3.Ds,obj.GD3.NumElt);

	% AFDiracDelta
	[h1x,h1y,h1z] = obj.GD3.Gradient(obj.AHeaviside);
	[h2x,h2y,h2z] = obj.GD3.Gradient(obj.Heaviside);
	[chx,chy,chz] = obj.GD3.CrossProduct(h2x,h2y,h2z,h1x,h1y,h1z);

	obj.AFDiracDelta = obj.GD3.DotProduct(obj.Cx,obj.Cy,obj.Cz,chx,chy,chz) ./ obj.NormCrossAF.^2;

	
end
