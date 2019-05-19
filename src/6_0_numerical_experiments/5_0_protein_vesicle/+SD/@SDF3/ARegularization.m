% calculate regularization terms for auxilary level set function
% rgNormal diffuses A field in the N direction
% rgSD makes A a signed distance function
% if UseDoubleWellPotential is true, then use a double well potential to
% regularize A field, otherwise use qudratic potential

function [rgNormal, rgSD] = ARegularization(obj, UseDoubleWellPotential)

	% calcualte gradient of Normal vectors
	tmp = obj.GD3.DotProduct(obj.Fxx, obj.Fxy, obj.Fzx, obj.Fx, obj.Fy, obj.Fz) ./ obj.FGradMag.^3;
	Nxx = obj.Fxx ./ obj.FGradMag - tmp .*obj.Fx;
	Nxy = obj.Fxy ./ obj.FGradMag - tmp .*obj.Fy;
	Nxz = obj.Fzx ./ obj.FGradMag - tmp .*obj.Fz;

	tmp = obj.GD3.DotProduct(obj.Fxy, obj.Fyy, obj.Fyz, obj.Fx, obj.Fy, obj.Fz) ./ obj.FGradMag.^3;
	Nyx = obj.Fxy ./ obj.FGradMag - tmp .*obj.Fx;
	Nyy = obj.Fyy ./ obj.FGradMag - tmp .*obj.Fy;
	Nyz = obj.Fyz ./ obj.FGradMag - tmp .*obj.Fz;

	tmp = obj.GD3.DotProduct(obj.Fzx, obj.Fyz, obj.Fzz, obj.Fx, obj.Fy, obj.Fz) ./ obj.FGradMag.^3;
	Nzx = obj.Fzx ./ obj.FGradMag - tmp .*obj.Fx;
	Nzy = obj.Fyz ./ obj.FGradMag - tmp .*obj.Fy;
	Nzz = obj.Fzz ./ obj.FGradMag - tmp .*obj.Fz;

	% calculate rgNormal
	rgNormal = obj.MeanCurvature .* obj.GD3.DotProduct(obj.Ax, obj.Ay, obj.Az, obj.Nx, obj.Ny, obj.Nz) ...
		+ obj.Nx .* obj.GD3.DotProduct(Nxx, Nyx, Nzx, obj.Ax, obj.Ay, obj.Az) ...
		+ obj.Ny .* obj.GD3.DotProduct(Nxy, Nyy, Nzy, obj.Ax, obj.Ay, obj.Az) ...
		+ obj.Nz .* obj.GD3.DotProduct(Nxz, Nyz, Nzz, obj.Ax, obj.Ay, obj.Az) ...
		+ obj.Nx .* obj.GD3.DotProduct(obj.Axx, obj.Axy, obj.Azx, obj.Nx, obj.Ny, obj.Nz) ...
		+ obj.Ny .* obj.GD3.DotProduct(obj.Axy, obj.Ayy, obj.Ayz, obj.Nx, obj.Ny, obj.Nz) ...
		+ obj.Nz .* obj.GD3.DotProduct(obj.Azx, obj.Ayz, obj.Azz, obj.Nx, obj.Ny, obj.Nz) ;

	% calculate rgSD
	pprime = obj.AGradMag - 1;
	pdprime = ones(obj.GD3.Size, 'gpuArray');
	if UseDoubleWellPotential
		Mask = obj.AGradMag > 1;
		pprime(Mask) = sin(2*pi*obj.AGradMag(Mask)) / (2*pi);
		pdprime(Mask) = cos(2*pi*obj.AGradMag(Mask));
	end
	rgSD = pprime .* obj.ALaplacian ./ obj.AGradMag ...
		+ (pdprime .* obj.AGradMag - pprime) ...
			.* ( obj.Ax .* obj.GD3.DotProduct(obj.Ax, obj.Ay, obj.Az, obj.Axx, obj.Axy, obj.Azx) + ...
				 obj.Ay .* obj.GD3.DotProduct(obj.Ax, obj.Ay, obj.Az, obj.Axy, obj.Ayy, obj.Ayz) + ...
				 obj.Az .* obj.GD3.DotProduct(obj.Ax, obj.Ay, obj.Az, obj.Azx, obj.Ayz, obj.Azz) ) ...
			./ obj.AGradMag.^3;



end
