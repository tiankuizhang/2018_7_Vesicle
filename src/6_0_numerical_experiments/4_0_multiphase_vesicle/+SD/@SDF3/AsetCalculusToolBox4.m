% calculate geometies of the Auxilary level set function, i.e. embedded curve geometries 
function AsetCalculusToolBox4(obj)

	[obj.Ax,obj.Ay,obj.Az] = obj.GD3.Gradient4(obj.A);
	obj.AGradMag = obj.GD3.Norm(obj.Ax,obj.Ay,obj.Az);

	[obj.Cx,obj.Cy,obj.Cz] = obj.GD3.CrossProduct(obj.Fx,obj.Fy,obj.Fz,obj.Ax,obj.Ay,obj.Az);
	obj.NormCrossAF = obj.GD3.Norm(obj.Cx,obj.Cy,obj.Cz);

	obj.tx = obj.Cx ./ obj.NormCrossAF;
	obj.ty = obj.Cy ./ obj.NormCrossAF;
	obj.tz = obj.Cz ./ obj.NormCrossAF;

	[obj.nx,obj.ny,obj.nz] = obj.GD3.CrossProduct(obj.tx,obj.ty,obj.tz,obj.Nx,obj.Ny,obj.Nz);

	[obj.Axx,obj.Ayy,obj.Azz,obj.Axy,obj.Ayz,obj.Azx] = obj.GD3.Hessian4(obj.A);
	obj.ALaplacian = obj.Axx + obj.Ayy + obj.Azz;

	% geodesic curvature
	vx = obj.tx.*obj.Fxx + obj.ty.*obj.Fxy + obj.tz.*obj.Fzx;
	vy = obj.tx.*obj.Fxy + obj.ty.*obj.Fyy + obj.tz.*obj.Fyz;
	vz = obj.tx.*obj.Fzx + obj.ty.*obj.Fyz + obj.tz.*obj.Fzz;

	[w1x,w1y,w1z] = obj.GD3.CrossProduct(vx,vy,vz,obj.Ax,obj.Ay,obj.Az);

	vx = obj.tx.*obj.Axx + obj.ty.*obj.Axy + obj.tz.*obj.Azx;
	vy = obj.tx.*obj.Axy + obj.ty.*obj.Ayy + obj.tz.*obj.Ayz;
	vz = obj.tx.*obj.Azx + obj.ty.*obj.Ayz + obj.tz.*obj.Azz;

	[w2x,w2y,w2z] = obj.GD3.CrossProduct(obj.Fx,obj.Fy,obj.Fz,vx,vy,vz);

	obj.GeodesicCurvature = obj.GD3.DotProduct(obj.nx,obj.ny,obj.nz,w1x+w2x,w1y+w2y,w1z+w2z) ./ obj.NormCrossAF;
	%obj.GeodesicCurvature = min(abs(obj.GeodesicCurvature), 5) .* sign(obj.GeodesicCurvature);

	%obj.NormalCurvature = obj.GD3.DotProduct(obj.Nx,obj.Ny,obj.Nz,w1x+w2x,w1y+w2y,w1z+w2z) ./ NormCrossAF;

	%% NormalCurvature,GeodesicTorsion, BPerpendicular
	% Gradient of Normal vectors: Nxy = partial_Nx/partial_y
	Nxx = obj.Fxx./obj.FGradMag - obj.Fx.*(obj.Fxx.*obj.Fx+obj.Fxy.*obj.Fy+obj.Fzx.*obj.Fz) ./ obj.FGradMag.^3;
	Nyx = obj.Fxy./obj.FGradMag - obj.Fy.*(obj.Fxx.*obj.Fx+obj.Fxy.*obj.Fy+obj.Fzx.*obj.Fz) ./ obj.FGradMag.^3;
	Nzx = obj.Fzx./obj.FGradMag - obj.Fz.*(obj.Fxx.*obj.Fx+obj.Fxy.*obj.Fy+obj.Fzx.*obj.Fz) ./ obj.FGradMag.^3;

	Nxy = obj.Fxy./obj.FGradMag - obj.Fx.*(obj.Fxy.*obj.Fx+obj.Fyy.*obj.Fy+obj.Fyz.*obj.Fz) ./ obj.FGradMag.^3;
	Nyy = obj.Fyy./obj.FGradMag - obj.Fy.*(obj.Fxy.*obj.Fx+obj.Fyy.*obj.Fy+obj.Fyz.*obj.Fz) ./ obj.FGradMag.^3;
	Nzy = obj.Fyz./obj.FGradMag - obj.Fz.*(obj.Fxy.*obj.Fx+obj.Fyy.*obj.Fy+obj.Fyz.*obj.Fz) ./ obj.FGradMag.^3;

	Nxz = obj.Fzx./obj.FGradMag - obj.Fx.*(obj.Fzx.*obj.Fx+obj.Fyz.*obj.Fy+obj.Fzz.*obj.Fz) ./ obj.FGradMag.^3;
	Nyz = obj.Fyz./obj.FGradMag - obj.Fy.*(obj.Fzx.*obj.Fx+obj.Fyz.*obj.Fy+obj.Fzz.*obj.Fz) ./ obj.FGradMag.^3;
	Nzz = obj.Fzz./obj.FGradMag - obj.Fz.*(obj.Fzx.*obj.Fx+obj.Fyz.*obj.Fy+obj.Fzz.*obj.Fz) ./ obj.FGradMag.^3;

	% NormalCurvature. (vx,vy,vz): gradient of Normal vectors along the curve
	vx =  Nxx .* obj.tx + Nxy .* obj.ty + Nxz .* obj.tz;
	vy =  Nyx .* obj.tx + Nyy .* obj.ty + Nyz .* obj.tz;
	vz =  Nzx .* obj.tx + Nzy .* obj.ty + Nzz .* obj.tz;

	obj.NormalCurvature = - obj.GD3.DotProduct(obj.tx,obj.ty,obj.tz,vx,vy,vz);

	% GeodesicTorsion, BPerpendicular. (vx,vy,vz): gradient of Normal vectors normal to the curve
	vx =  Nxx .* obj.nx + Nxy .* obj.ny + Nxz .* obj.nz;
	vy =  Nyx .* obj.nx + Nyy .* obj.ny + Nyz .* obj.nz;
	vz =  Nzx .* obj.nx + Nzy .* obj.ny + Nzz .* obj.nz;
							  
	obj.GeodesicTorsion = - obj.GD3.DotProduct(obj.tx,obj.ty,obj.tz,vx,vy,vz);
	obj.BPerpendicular  = - obj.GD3.DotProduct(obj.nx,obj.ny,obj.nz,vx,vy,vz);

	%%  primal of Heaviside(A), Heaviside(A), DiracDelta(A)
	obj.AHPrimal = max(obj.A,0);

	[Px,Py,Pz] = obj.GD3.Gradient4(obj.AHPrimal);
	P_lap = obj.GD3.Laplacian4(obj.AHPrimal);

	dot_DAHPrimal_DA = obj.GD3.DotProduct(Px,Py,Pz,obj.Ax,obj.Ay,obj.Az);

	obj.AHeaviside = dot_DAHPrimal_DA ./ obj.AGradMag.^2;
	obj.ADiracDelta = P_lap ./ obj.AGradMag.^2 - ...
		dot_DAHPrimal_DA .* obj.ALaplacian ./ obj.AGradMag.^4;

	% AFDiracDelta
	[h1x,h1y,h1z] = obj.GD3.Gradient4(obj.AHeaviside);
	[h2x,h2y,h2z] = obj.GD3.Gradient4(obj.Heaviside);
	[chx,chy,chz] = obj.GD3.CrossProduct(h2x,h2y,h2z,h1x,h1y,h1z);
	obj.AFDiracDelta = obj.GD3.DotProduct(obj.Cx,obj.Cy,obj.Cz,chx,chy,chz) ./ obj.NormCrossAF.^2;

	obj.ADiracDeltaDn = (obj.GD3.Laplacian4(obj.AHeaviside) - obj.ADiracDelta .* obj.GD3.Laplacian4(obj.A) ) ./ obj.AGradMag;
end









