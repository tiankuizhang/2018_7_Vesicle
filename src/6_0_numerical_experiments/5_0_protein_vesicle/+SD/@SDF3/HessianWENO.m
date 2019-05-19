function [fxx,fyy,fzz,fxy,fyz,fzx] = HessianWENO(obj,f)

	[fx,fy,fz] = obj.GradientWENO(f);

	[fxx,fxy,fxz] = obj.GradientWENO(fx);
	[fyx,fyy,fyz] = obj.GradientWENO(fy);
	[fzx,fzy,fzz] = obj.GradientWENO(fz);

	fxy = (fxy + fyx) / 2.0;
	fyz = (fyz + fzy) / 2.0;
	fzx = (fzx + fxz) / 2.0;
end
