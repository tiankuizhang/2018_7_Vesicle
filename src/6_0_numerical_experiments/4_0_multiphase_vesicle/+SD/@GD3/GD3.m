classdef GD3 < handle
	
	%Grid3 : meshgrid in 3D

	properties (SetAccess = immutable)
		X % meshgrid X
		Y % meshgrid Y
		Z % meshgrid Z
		Dx % spacing in x direction
		Dy % spacing in y direction
		Dz % spacing in z direction
		Ds % mean spacing
		ncols % number of grid points in x direction
		mrows % number of grid points in y direction
		lshts % number of grid points in z direction
		Size % size array
		NumElt % total number of grid points 
		xmin % min value of coordinate x 
		xmax % max value of coordinate x
		ymin % min value of coordinate y
		ymax % max value of coordinate y
		zmin % min value of coordinate z
		zmax % max value of coordinate z
		BOX % bounding box of grid domain
		
		% index matrix and the shifted index matrices
		ooo % not shifted
		
		oXo % shifted in the -x direction thus represent indices to the right  
		oxo % shifted in the +x direction thus represent indices to the left
		Yoo 
		yoo
		ooZ
		ooz

		YXo
		yXo
		Yxo
		yxo
		
		YoZ
		Yoz
		yoZ
		yoz
		
		oXZ
		oxZ
		oXz
		oxz
		
		YXZ
		YXz
		YxZ
		Yxz
		yXZ
		yXz
		yxZ
		yxz

		% first/sencond order derivative operators (including cross derivatives)
		% using centeral difference with periodic boundary condition
		Lx
		Ly
		Lz
		Lxx
		Lyy
		Lzz
		Lxy
		Lyz
		Lzx
		LLaplacian
		LBiLaplacian

		% identity matrix
		Idt

		% frequency domain
		kx
		ky
		kz

	end

	methods

		function obj = GD3(Xm, Ym, Zm)
			obj.X = Xm;
			obj.Y = Ym;
			obj.Z = Zm;
			obj.Dx = Xm(1,2,1) - Xm(1,1,1);
			obj.Dy = Ym(2,1,1) - Ym(1,1,1);
			obj.Dz = Zm(1,1,2) - Zm(1,1,1);
			obj.Ds = (obj.Dx + obj.Dy + obj.Dz) / 3;
			obj.Size = size(Xm);
			obj.mrows = obj.Size(1);
			obj.ncols = obj.Size(2);
			obj.lshts = obj.Size(3);
			obj.NumElt = prod(obj.Size);
			obj.xmin = Xm(1,1,1);
			obj.ymin = Ym(1,1,1);
			obj.zmin = Zm(1,1,1);
			obj.xmax = Xm(1,end,1);
			obj.ymax = Ym(end,1,1);
			obj.zmax = Zm(1,1,end);
			obj.BOX = gather( [obj.xmin obj.xmax obj.ymin obj.ymax obj.zmin obj.zmax] );

			% the largest integer uint32 can represent is 2^32~4e9~(1.6e3)^3
			% the largest integer uint16 can represent is 2^16~6e4~(4e2)^40^3
			% thus as long as NumElt<2^32 we can use uint32 
			obj.ooo = gpuArray(reshape(1:obj.NumElt,obj.Size));

			obj.oXo = circshift(obj.ooo, [	0	-1	0	]);
			obj.oxo = circshift(obj.ooo, [	0	1	0	]);
			obj.Yoo = circshift(obj.ooo, [	-1	0	0	]);
			obj.yoo = circshift(obj.ooo, [	1	0	0	]);
			obj.ooZ = circshift(obj.ooo, [	0	0	-1	]);
			obj.ooz = circshift(obj.ooo, [	0	0	1	]);

			obj.YXo = circshift(obj.ooo, [	-1	-1	0	]);
			obj.yXo = circshift(obj.ooo, [	1	-1	0	]); 
			obj.Yxo = circshift(obj.ooo, [	-1	1	0	]);
			obj.yxo = circshift(obj.ooo, [	1	1	0	]); 

			obj.YoZ = circshift(obj.ooo, [	-1	0	-1	]);
			obj.Yoz = circshift(obj.ooo, [	-1	0	1	]); 
			obj.yoZ = circshift(obj.ooo, [	1	0	-1	]);
			obj.yoz = circshift(obj.ooo, [	1	0	1	]); 

			obj.oXZ = circshift(obj.ooo, [	0	-1	-1	]);
			obj.oxZ = circshift(obj.ooo, [	0	1	-1	]); 
			obj.oXz = circshift(obj.ooo, [	0	-1	1	]);
			obj.oxz = circshift(obj.ooo, [	0	1	1	]); 

			obj.YXZ = circshift(obj.ooo, [	-1	-1	-1	]);
			obj.YXz = circshift(obj.ooo, [	-1	-1	1	]); 
			obj.YxZ = circshift(obj.ooo, [	-1	1	-1	]);
			obj.Yxz = circshift(obj.ooo, [	-1	1	1	]); 
			obj.yXZ = circshift(obj.ooo, [	1	-1	-1	]);
			obj.yXz = circshift(obj.ooo, [	1	-1	1	]); 
			obj.yxZ = circshift(obj.ooo, [	1	1	-1	]);
			obj.yxz = circshift(obj.ooo, [	1	1	1	]); 

			obj.Lx  =	( sparse(obj.ooo(:), obj.oXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxo(:),	-1,	obj.NumElt, obj.NumElt) ) ...
		   				* (1 / (2*obj.Dx));

			obj.Ly  =	( sparse(obj.ooo(:), obj.Yoo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoo(:),	-1,	obj.NumElt, obj.NumElt) ) ...
		  			 	* (1 / (2*obj.Dy));

			obj.Lz  = 	( sparse(obj.ooo(:), obj.ooZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooz(:),	-1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (2*obj.Dz)); 

			obj.Lxx =	( sparse(obj.ooo(:), obj.oXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxo(:),	1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (obj.Dx.^2));

			obj.Lyy =	( sparse(obj.ooo(:), obj.Yoo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoo(:),	1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (obj.Dy.^2));

			obj.Lzz = 	( sparse(obj.ooo(:), obj.ooZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooo(:),	-2,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.ooz(:),	1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (obj.Dz.^2));

			obj.Lxy = 	( sparse(obj.ooo(:), obj.YXo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yxo(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.Yxo(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yXo(:),	-1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (4*obj.Ds.^2)); 

			obj.Lyz = 	( sparse(obj.ooo(:), obj.YoZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoz(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.Yoz(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.yoZ(:),	-1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (4*obj.Ds.^2)); 

			obj.Lzx = 	( sparse(obj.ooo(:), obj.oXZ(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxz(:),	1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oXz(:),	-1,	obj.NumElt, obj.NumElt) ...
						+ sparse(obj.ooo(:), obj.oxZ(:),	-1,	obj.NumElt, obj.NumElt) ) ...
						* (1 / (4*obj.Ds.^2)); 

			obj.LLaplacian = obj.Lxx + obj.Lyy + obj.Lzz;
			obj.LBiLaplacian = obj.LLaplacian * obj.LLaplacian;

			obj.Idt = sparse(obj.ooo(:), obj.ooo(:), 1, obj.NumElt, obj.NumElt);

			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% set coordinate for the frequency domain
			kxcpu = 2 * pi * (-obj.ncols/2 : obj.ncols/2-1) / (obj.ncols * obj.Dx);
			kycpu = 2 * pi * (-obj.mrows/2 : obj.mrows/2-1) / (obj.mrows * obj.Dy);
			kzcpu = 2 * pi * (-obj.lshts/2 : obj.lshts/2-1) / (obj.lshts * obj.Dz);

			kxcpu = fftshift(kxcpu);
			kycpu = fftshift(kycpu);
			kzcpu = fftshift(kzcpu);

			[kx, ky, kz] = meshgrid(kxcpu, kycpu, kzcpu);

			obj.kx = gpuArray(kx);
			obj.ky = gpuArray(ky);
			obj.kz = gpuArray(kz);
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		end

		% create a sparse diagonal matrix out of a field
		% look up the function spdiags and replace it with the built-in function
		function val = SparseDiag(obj, Field)
			val = sparse(obj.ooo(:), obj.ooo(:), reshape(Field,[],1), obj.NumElt, obj.NumElt);
		end

		% various derivatives of a field
		function val = Fx(obj, Field)
			val = (Field(obj.oXo) - Field(obj.oxo)) / (2*obj.Dx);
		end

		function val = Fx4(obj, Field)
			F2  = circshift(Field,[0,-2,0]); 
			Fm2 = circshift(Field,[0, 2,0]);
			val = (-F2 + 8*Field(obj.oXo) - 8*Field(obj.oxo) + Fm2) / (12*obj.Dx);
		end

		function val = Fy(obj, Field)
			val = (Field(obj.Yoo) - Field(obj.yoo)) / (2*obj.Dy);
		end

		function val = Fy4(obj, Field)
			F2  = circshift(Field,[-2,0,0]);
			Fm2 = circshift(Field,[ 2,0,0]);
			val = (-F2 + 8*Field(obj.Yoo) - 8*Field(obj.yoo) + Fm2) / (12*obj.Dy); 
		end

		function val = Fz(obj, Field)
			val = (Field(obj.ooZ) - Field(obj.ooz)) / (2*obj.Dz);
		end

		function val = Fz4(obj, Field)
			F2  = circshift(Field,[0,0,-2]);
			Fm2 = circshift(Field,[0,0, 2]);
			val = (-F2 + 8*Field(obj.ooZ) - 8*Field(obj.ooz) + Fm2) / (12*obj.Dz); 
		end

		function val = Fxx(obj, Field)
			val = (Field(obj.oXo) - 2*Field + Field(obj.oxo)) / (obj.Dx.^2);
		end

		function val = Fxx4(obj, Field)
			F2  = circshift(Field,[0,-2,0]); 
			Fm2 = circshift(Field,[0, 2,0]);
			val = (-F2 + 16*Field(obj.oXo) -30*Field + 16*Field(obj.oxo) - Fm2) / (12*obj.Dx.^2);
		end

		function val = Fyy(obj, Field)
			val = (Field(obj.Yoo) - 2*Field + Field(obj.yoo)) / (obj.Dy.^2);
		end

		function val = Fyy4(obj, Field)
			F2  = circshift(Field,[-2,0,0]);
			Fm2 = circshift(Field,[ 2,0,0]);
			val = (-F2 + 16*Field(obj.Yoo) - 30*Field + 16*Field(obj.yoo) - Fm2) / (12*obj.Dy.^2); 
		end

		function val = Fzz(obj, Field)
			val = (Field(obj.ooZ) - 2*Field + Field(obj.ooz)) / (obj.Dz.^2);
		end

		function val = Fzz4(obj, Field)
			F2  = circshift(Field,[0,0,-2]);
			Fm2 = circshift(Field,[0,0, 2]);
			val = (-F2 + 16*Field(obj.ooZ) - 30*Field + 16*Field(obj.ooz) - Fm2) / (12*obj.Dz.^2); 
		end

		function val = Fxy(obj, Field)
			val = (Field(obj.YXo) + Field(obj.yxo) - Field(obj.Yxo) - Field(obj.yXo)) / (4*obj.Ds.^2);
		end

		function val = Fxy4(obj, Field)
			F22 = circshift(Field,[-2,-2,0]);
			Fm2m2 = circshift(Field,[2,2,0]);
			F2m2 = circshift(Field,[-2,2,0]);
			Fm22 = circshift(Field,[2,-2,0]);
			val = (-F22 - Fm2m2 + F2m2 + Fm22 + 16*Field(obj.YXo) + 16 * Field(obj.yxo) - 16 * Field(obj.Yxo) - 16 * Field(obj.yXo)) / (48*obj.Dx*obj.Dy);
		end	

		function val = Fyz(obj, Field)
			val = (Field(obj.YoZ) + Field(obj.yoz) - Field(obj.Yoz) - Field(obj.yoZ)) / (4*obj.Ds.^2);
		end

		function val = Fyz4(obj, Field)
			F22 = circshift(Field,[-2,0,-2]);
			Fm2m2 = circshift(Field,[2,0,2]);
			F2m2 = circshift(Field,[-2,0,2]);
			Fm22 = circshift(Field,[2,0,-2]);
			val = (-F22 - Fm2m2 + F2m2 + Fm22 + 16*Field(obj.YoZ) + 16 * Field(obj.yoz) - 16 * Field(obj.Yoz) - 16 * Field(obj.yoZ)) / (48*obj.Dz*obj.Dy);
		end	

		function val = Fzx(obj, Field)
			val = (Field(obj.oXZ) + Field(obj.oxz) - Field(obj.oXz) - Field(obj.oxZ)) / (4*obj.Ds.^2);
		end

		function val = Fzx4(obj, Field)
			F22 = circshift(Field,[0,-2,-2]);
			Fm2m2 = circshift(Field,[0,2,2]);
			F2m2 = circshift(Field,[0,2,-2]);
			Fm22 = circshift(Field,[0,-2,2]);
			val = (-F22 - Fm2m2 + F2m2 + Fm22 + 16*Field(obj.oXZ) + 16 * Field(obj.oxz) - 16 * Field(obj.oXz) - 16 * Field(obj.oxZ)) / (48*obj.Dz*obj.Dx);
		end	

		function [fx,fy,fz] = Gradient(obj,f)
			fx = obj.Fx(f);
			fy = obj.Fy(f);
			fz = obj.Fz(f);
		end

		function [fx,fy,fz] = Gradient4(obj,f)
			fx = obj.Fx4(f);
			fy = obj.Fy4(f);
			fz = obj.Fz4(f);
		end

		function [fxx,fyy,fzz,fxy,fyz,fzx] = Hessian(obj,f)
			fxx = obj.Fxx(f);
			fyy = obj.Fyy(f);
			fzz = obj.Fzz(f);
			fxy = obj.Fxy(f);
			fyz = obj.Fyz(f);
			fzx = obj.Fzx(f);
		end

		function [fxx,fyy,fzz,fxy,fyz,fzx] = Hessian4(obj,f)
			fxx = obj.Fxx4(f);
			fyy = obj.Fyy4(f);
			fzz = obj.Fzz4(f);
			fxy = obj.Fxy4(f);
			fyz = obj.Fyz4(f);
			fzx = obj.Fzx4(f);
		end

		function val = Laplacian(obj, Field)
			val = obj.Fxx(Field) + obj.Fyy(Field) + obj.Fzz(Field);
		end

		function val = Laplacian4(obj, Field)
			val = obj.Fxx4(Field) + obj.Fyy4(Field) + obj.Fzz4(Field);
		end

		function [wx,wy,wz] = CrossProduct(obj,ux,uy,uz,vx,vy,vz)
			wx = uy.*vz - uz.*vy;
			wy = uz.*vx - ux.*vz;
			wz = ux.*vy - uy.*vx;
		end

		function val = DotProduct(obj,ux,uy,uz,vx,vy,vz)
			val = ux.*vx + uy.*vy + uz.*vz;
		end

		function mag = Norm(obj,vx,vy,vz)
			mag = max(sqrt(vx.^2+vy.^2+vz.^2),1e-14);	
		end
		
	end

	methods
		[WENO_back_x, WENO_fore_x, WENO_back_y, WENO_fore_y, WENO_back_z, WENO_fore_z] = WENODerivative(obj, Field);
	end

	properties
		% gpu Kernel function for the calculation of weno derivatives
		weno_derivative
	end
	methods
		[WENO_back_x, WENO_fore_x, WENO_back_y, WENO_fore_y, WENO_back_z, WENO_fore_z] = GPUWENODerivative(obj, Field);
	end

	methods
		function DrawBox2(obj)
			hold on
			line([obj.xmin,obj.xmax,obj.xmax,obj.xmin,obj.xmin],...
				 [obj.ymin,obj.ymin,obj.ymax,obj.ymax,obj.ymin],...
				 [obj.zmin,obj.zmin,obj.zmin,obj.zmin,obj.zmin],...
				 'Color','k');
			line([obj.xmin,obj.xmax,obj.xmax,obj.xmin,obj.xmin],...
				 [obj.ymin,obj.ymin,obj.ymax,obj.ymax,obj.ymin],...
				 [obj.zmax,obj.zmax,obj.zmax,obj.zmax,obj.zmax], ...
				 'Color','k');
			line([obj.xmin,obj.xmin],[obj.ymin,obj.ymin],[obj.zmin,obj.zmax],'Color','k')
			line([obj.xmin,obj.xmin],[obj.ymax,obj.ymax],[obj.zmin,obj.zmax],'Color','k')
			line([obj.xmax,obj.xmax],[obj.ymin,obj.ymin],[obj.zmin,obj.zmax],'Color','k')
			line([obj.xmax,obj.xmax],[obj.ymax,obj.ymax],[obj.zmin,obj.zmax],'Color','k')
			hold off
		end
		function DrawBox(obj)
			hold on
			%line([obj.xmin,obj.xmax,obj.xmax,obj.xmax,obj.xmin,obj.xmin,obj.xmin],...
			%	 [obj.ymin,obj.ymin,obj.ymin,obj.ymax,obj.ymax,obj.ymax,obj.ymin],...
			%	 [obj.zmin,obj.zmin,obj.zmax,obj.zmax,obj.zmax,obj.zmin,obj.zmin],...
			%	 'Color','k');
			line([obj.xmin,obj.xmax,obj.xmax,obj.xmax,obj.xmin,obj.xmin,obj.xmin],...
				 [obj.ymin,obj.ymin,obj.ymax,obj.ymax,obj.ymax,obj.ymin,obj.ymin],...
				 [obj.zmin,obj.zmin,obj.zmin,obj.zmax,obj.zmax,obj.zmax,obj.zmin],...
				 'Color','k');
			hold off
		end
	end

end
