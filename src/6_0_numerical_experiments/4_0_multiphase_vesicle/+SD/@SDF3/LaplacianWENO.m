function val = LaplacianWENO(obj,f)

	[fxx,fyy,fzz,fxy,fyz,fzx] = obj.HessianWENO(f);

	val = fxx + fyy + fzz;
end
