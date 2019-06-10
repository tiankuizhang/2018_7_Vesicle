% calculate regularization terms for the level set function
% rgSD makes F a signed distance function
% if UseDoubleWellPotential is true, then use a double well potential to
% regularize A field, otherwise use qudratic potential

function rgSD = ARegularization(obj, UseDoubleWellPotential)

	% calculate rgSD
	pprime = obj.FGradMag - 1;
	pdprime = ones(obj.GD3.Size, 'gpuArray');
	if UseDoubleWellPotential
		Mask = obj.FGradMag > 1;
		pprime(Mask) = sin(2*pi*obj.FGradMag(Mask)) / (2*pi);
		pdprime(Mask) = cos(2*pi*obj.FGradMag(Mask));
	end
	rgSD = pprime .* obj.FLaplacian ./ obj.FGradMag ...
		+ (pdprime .* obj.FGradMag - pprime) ...
			.* ( obj.Fx .* obj.GD3.DotProduct(obj.Fx, obj.Fy, obj.Fz, obj.Fxx, obj.xy, obj.Fzx) + ...
				 obj.Fy .* obj.GD3.DotProduct(obj.Fx, obj.Fy, obj.Fz, obj.Fxy, obj.Fyy, obj.Fyz) + ...
				 obj.Fz .* obj.GD3.DotProduct(obj.Fx, obj.Fy, obj.Fz, obj.Fzx, obj.Fyz, obj.Fzz) ) ...
			./ obj.FGradMag.^3;



end
