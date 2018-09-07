% calculate gradient, normal, second derivatives, curvatures, Dirac_Delta function
% Heaviside function
function setCalculusToolBox(obj)

	[obj.Fx,obj.Fy,obj.Fz] = obj.GD3.Gradient(obj.F);
	obj.FGradMag = sqrt(obj.Fx.^2+obj.Fy.^2+obj.Fz.^2);

	obj.Nx = obj.Fx ./ obj.FGradMag;
	obj.Ny = obj.Fy ./ obj.FGradMag;
	obj.Nz = obj.Fz ./ obj.FGradMag;

	[obj.Fxx,obj.Fyy,obj.Fzz,obj.Fxy,obj.Fyz,obj.Fzx] = obj.GD3.Hessian(obj.F);
	obj.FLaplacian = obj.GD3.Laplacian(obj.F);

	% calculate mean curvature
	col1 = obj.Fxx.*obj.Fx + obj.Fxy.*obj.Fy + obj.Fzx.*obj.Fz;
	col2 = obj.Fxy.*obj.Fx + obj.Fyy.*obj.Fy + obj.Fyz.*obj.Fz;
	col3 = obj.Fzx.*obj.Fx + obj.Fyz.*obj.Fy + obj.Fzz.*obj.Fz;

	obj.MeanCurvature = obj.FLaplacian ./ obj.FGradMag - ...
		(obj.Fx.*col1 + obj.Fy.*col2 + obj.Fz.*col3) ./ obj.FGradMag.^3;
	
	% calculate the Gaussian Curvature
	col1 = (obj.Fyy.*obj.Fzz-obj.Fyz.*obj.Fyz).*obj.Fx + ...
		   (obj.Fzx.*obj.Fyz-obj.Fxy.*obj.Fzz).*obj.Fy + ...
		   (obj.Fxy.*obj.Fyz-obj.Fzx.*obj.Fyy).*obj.Fz ;
	col2 = (obj.Fyz.*obj.Fzx-obj.Fxy.*obj.Fzz).*obj.Fz + ...
		   (obj.Fxx.*obj.Fzz-obj.Fzx.*obj.Fzx).*obj.Fy + ...
		   (obj.Fzx.*obj.Fxy-obj.Fxx.*obj.Fyz).*obj.Fz ;
	col3 = (obj.Fxy.*obj.Fyz-obj.Fyy.*obj.Fzx).*obj.Fx + ...
		   (obj.Fzx.*obj.Fxy-obj.Fxx.*obj.Fyz).*obj.Fy + ...
		   (obj.Fxx.*obj.Fyy-obj.Fxy.*obj.Fxy).*obj.Fz ;

	obj.GaussianCurvature = ...
		(obj.Fx.*col1 + obj.Fy.*col2 + obj.Fz.*col3) ./ obj.FGradMag.^4;

	% primal of Heaviside, Heaviside, Dirac_Delta
	obj.HPrimal = max(obj.F,0);

	[Px,Py,Pz] = obj.GD3.Gradient(obj.HPrimal);

	P_lap = obj.GD3.Laplacian(obj.HPrimal);

	% dot product of gradient of HPrimal and F
	dot_DHPrimal_DF = obj.GD3.DotProduct(Px,Py,Pz,obj.Fx,obj.Fy,obj.Fz);

	obj.Heaviside = dot_DHPrimal_DF ./ obj.FGradMag.^2;

	obj.DiracDelta = P_lap ./ obj.FGradMag.^2 - ...
		dot_DHPrimal_DF .* obj.FLaplacian ./ obj.FGradMag.^4;




end
