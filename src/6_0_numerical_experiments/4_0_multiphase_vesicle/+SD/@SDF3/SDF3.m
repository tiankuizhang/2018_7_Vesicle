classdef SDF3 < handle
	
	%SDF3 : signed distance function in 3D

	properties (SetAccess = immutable)
		GD3 % SD.GD3 object
	end

	methods

		function obj = SDF3(grid, Xm, Ym, Zm, Val)
			obj.GD3 = grid;
			obj.GPUInitialize;
			obj.F = Val;
		end

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculus tool box : derivatives, curvature, DiracDelta function, Heaviside function
% 53 lines
	properties
		F % values of the signed distance function

		Fx % gradient of F
		Fy
		Fz
		FGradMag % magnitude of (Fx,Fy,Fz)

		Nx % normal of the level set function
		Ny
		Nz

		Fxx % second derivatives
		Fyy
		Fzz
		Fxy
		Fyz
		Fzx
		FLaplacian

		MeanCurvature
		GaussianCurvature

		HPrimal % max(0,F)
		Heaviside
		DiracDelta

		xpr % distance to the right boundary if there is boundary to the right
		xpl
		ypf
		ypb
		zpu
		zpd
		min_dist

		Box
	end


	methods
		% update the above properties
		setCalculusToolBox(obj)
		setCalculusToolBox4(obj)
		setCalculusToolBoxGA(obj,threshold)
		setCalculustToolBoxWENO(obj)
		GPUsetCalculusToolBox(obj)  
		setDistance(obj) % used after F modified and before any use of xpr etc.
	end

	methods
		% return area and volume of the closed surface
		function area = calArea(obj)
			mask = abs(obj.F) < 2*obj.GD3.Ds;
			area = sum(obj.DiracDelta(mask)) * obj.GD3.Ds.^3;
		end
		function volume = calVolume(obj)
			mask = obj.F < 2*obj.GD3.Ds;
			volume = sum(1-obj.Heaviside(mask)) * obj.GD3.Ds.^3;
		end
		% surface integral
		function integral = surfaceIntegral(obj,field)
			mask = abs(obj.F) < 2*obj.GD3.Ds;
			integral = sum(field(mask).*obj.DiracDelta(mask).*obj.FGradMag(mask))*obj.GD3.Ds.^3;
		end
		% surface Laplacian of field
		function val = SurfaceLaplacian(obj, field)
			[fxx,fyy,fzz,fxy,fyz,fzx] = obj.GD3.Hessian(field);
			val = fxx + fyy + fzz - ...
					- obj.Nx.^2 .* fxx ...
					- obj.Ny.^2 .* fyy ...
					- obj.Nz.^2 .* fzz ...
					- 2 * obj.Nx .* obj.Ny .* fxy ...
					- 2 * obj.Ny .* obj.Nz .* fyz ...
					- 2 * obj.Nz .* obj.Nx .* fzx ;
		end
		function val = SurfaceLaplacian4(obj, field)
			[fxx,fyy,fzz,fxy,fyz,fzx] = obj.GD3.Hessian4(field);
			val = fxx + fyy + fzz - ...
					- obj.Nx.^2 .* fxx ...
					- obj.Ny.^2 .* fyy ...
					- obj.Nz.^2 .* fzz ...
					- 2 * obj.Nx .* obj.Ny .* fxy ...
					- 2 * obj.Ny .* obj.Nz .* fyz ...
					- 2 * obj.Nz .* obj.Nx .* fzx ;
		end
		% surface divergence of vector field
		function val = SurfaceDivergence(obj, vx, vy, vz)
			[vxx, vxy, vxz] = obj.GD3.Gradient(vx);
			[vyx, vyy, vyz] = obj.GD3.Gradient(vy);
			[vzx, vzy, vzz] = obj.GD3.Gradient(vz);
			
			val = vxx + vyy + vzz ...
				- obj.Nx .* obj.GD3.DotProduct(obj.Nx, obj.Ny, obj.Nz, vxx, vxy, vxz) ...
				- obj.Ny .* obj.GD3.DotProduct(obj.Nx, obj.Ny, obj.Nz, vyx, vyy, vyz) ...
				- obj.Nz .* obj.GD3.DotProduct(obj.Nx, obj.Ny, obj.Nz, vzx, vzy, vzz);
		end

		% calculate local Lagrange multiplier constraining local area
		function [localTension,residual] = localLagrangeMultiplier(obj,rhs,Dt,Alpha)
			% operator to be solved
			Op = obj.GD3.Lxx + obj.GD3.Lyy + obj.GD3.Lzz - ...
					(   obj.GD3.SparseDiag(obj.Nx .^2) * obj.GD3.Lxx + ...
						obj.GD3.SparseDiag(obj.Ny .^2) * obj.GD3.Lyy + ...
						obj.GD3.SparseDiag(obj.Nz .^2) * obj.GD3.Lzz + ...
						obj.GD3.SparseDiag(obj.Nx .* obj.Ny) * obj.GD3.Lxy * 2 + ...
						obj.GD3.SparseDiag(obj.Ny .* obj.Nz) * obj.GD3.Lyz * 2 + ...
						obj.GD3.SparseDiag(obj.Nz .* obj.Nx) * obj.GD3.Lzx * 2 ...
					) ...
				- obj.GD3.SparseDiag(obj.MeanCurvature.^2) ;
			% right hand side
			S = reshape(rhs, [obj.GD3.NumElt, 1]);
			% precompute division
			KSquaredInverse = -1./(obj.GD3.kx.^2 + obj.GD3.ky.^2 + obj.GD3.kz.^2 + Alpha);
			% it seems that restart about {10,...,15} gives best speed
			%tic
			%[localTension,~,~,~,~] = gmres(Op, S, 11, 1e-6, 300, @mfun);
			%localTension = gmres(Op, S, 11, 1e-6, 300, @mfun);
			%localTension = gmres(Op, S, 20, 1e-1, 300, @mfun);
			[localTension,~,~,~,~] = gmres(Op, S, 11, 1e-2, 300, @mfun);
			%toc
			residual = Op * localTension - S;
			localTension = reshape(localTension, obj.GD3.Size);
			residual = reshape(residual, obj.GD3.Size);
			% preconditioner, 
			% Alpha = (c22+c33+c23+c32)/c11 as the mean squared MeanCurvature
			function y = mfun(S)
				fftS = fftn(reshape(S,obj.GD3.Size)) .* KSquaredInverse;
				y = reshape(real(ifftn(fftS)), [obj.GD3.NumElt, 1]);
			end
		end

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	properties
		A % values of the auxilary level set function
		
		Ax
		Ay
		Az
		AGradMag

		Cx % cross product of grad(F) and grad(A)
		Cy
		Cz
		NormCrossAF % norm of C

		tx % tangent vectors of the curve
		ty
		tz

		nx % surface normal of the curve
		ny
		nz

		Axx % Hessian of A
		Ayy
		Azz
		Axy
		Ayz
		Azx
		ALaplacian

		GeodesicCurvature
		NormalCurvature
		GeodesicTorsion
		BPerpendicular

		AHPrimal
		AHeaviside
		ADiracDelta
		ADiracDeltaDn
		AFDiracDelta % product of ADiracDelta and FDiracDelta

	end

	methods
		AsetCalculusToolBox(obj)
		AsetCalculusToolBox4(obj)
		GPUAsetCalculusToolBox(obj)
	end

	methods
		% return length of curve
		function length = calLength(obj)
			mask = abs(obj.F) < 2*obj.GD3.Ds & abs(obj.A) < 2*obj.GD3.Ds;
			length = sum(obj.AFDiracDelta(mask).*obj.NormCrossAF(mask)) * obj.GD3.Ds.^3;
		end
		function val = LineIntegral(obj, field)
			mask = abs(obj.F) < 2*obj.GD3.Ds & abs(obj.A) < 2*obj.GD3.Ds;
			val = sum(field(mask).*obj.AFDiracDelta(mask).*obj.NormCrossAF(mask)) * obj.GD3.Ds.^3;
		end
		function area = AcalArea(obj)
			mask = abs(obj.F) < 2*obj.GD3.Ds;
			area = sum( obj.DiracDelta(mask).*(1-obj.AHeaviside(mask)) ) * obj.GD3.Ds.^3;
		end
		function integral = AsurfaceIntegral(obj,field)
			mask = abs(obj.F) < 2*obj.GD3.Ds;
			integral = sum(field(mask).*obj.DiracDelta(mask).*obj.FGradMag(mask).*(1-obj.AHeaviside(mask)))*obj.GD3.Ds.^3;
		end
		function val = BivalueField(obj, valNeg, valPos)
			val = valPos * obj.AHeaviside + valNeg * (1 - obj.AHeaviside) ;
			%val = valNeg * ones(obj.GD3.Size,'gpuArray');
			%val(obj.A>0) = valPos;
		end
		function [ft,fn] = AsurfaceDerivative(obj,field)
			% assuming field is alread extended away from the surface
			[fx,fy,fz] = obj.GD3.Gradient(field);
			ft = obj.GD3.DotProduct(fx,fy,fz,obj.tx,obj.ty,obj.tz);
			fn = obj.GD3.DotProduct(fx,fy,fz,obj.nx,obj.ny,obj.nz);
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GPU related properties and functions: 40 lines
	properties
		
		% parameter for GPU kernel
			ThreadBlockSize
			GridSize

		% kernel funcions object for ENORK2 reinitialization scheme
			% calculate grid step with modification near boundary
			ENORK2_boundary_correction 
			% calculate the numerical Hamiltonian for the Reinitalization equation
			ENORK2_reinitiliaztion_step 

		% kernel function object for WENORK3 reinitialization  scheme
			% calculate distance to the nearby nodes including boundary nodes
			cubic_boundary_correction
			% calculate the reinitialization step
			WENORK3_reinitialization_step

		% kernel function object WENO5RK3 reinitialization scheme
			% calculate distance to the bounary along x,y,z directions
			sixth_boundary_correction
			% calculate the reinitialization step
			WENO5RK3_reinitialization_step

		% kernel functions object for ENORK2 extend scheme
			ENORK2_upwind_normal % calculate upwind normals of the level set function
			ENORK2_boundary_interpolate % interpolate values at the boundary
			ENORK2_extend_step % calculate the extension step

		% kernel functions object for WENORK3 extend scheme
			WENORK3_upwind_normal % calculate upwind normals of the level set function
			WENORK3_boundary_interpolate % interpolate values at the boundary
			WENORK3_extend_step % calculate the extension step

		% kernel functions object for WENORK3 extend scheme
			WENO5RK3_upwind_normal % calculate upwind normals of the level set function
			WENO5RK3_boundary_interpolant % interpolate values at the boundary
			WENO5RK3_extend_step % calculate the extension step

		% kernel function object for GPUsetCalculusToolBox scheme
			set_calculus_toolbox % set Fx,Fy ...
			auxi_set_calculus_toolbox % set Ax,Ay ...
			upwind_derivative % calculate gradient of some field
			ga_set_calculus_toolbox % geometry aware version of set_calculus_toolbox
			advection_step % calculate advection term upwindly with WENO scheme

		% kernel function for numerical Hamiltonian for surface consevation law
			surface_conservation_step
			spatial_finite_volume_step

	end
	
	methods 
		GPUInitialize(obj)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilities : reinitliazation, extend, surface_redistance
	methods
		% reinitialization schemes of various orders of accuracy
		NewF = ENORK2Reinitialization(obj,F,iteration)
		NewF = WENORK3Reinitialization(obj,F,iteration)
		NewF = WENO5RK3Reinitialization(obj,F,iteration)
		
		% extend schems of various orders of accuracy
		NewC = ENORK2Extend(obj, C, iteration)
		NewC = WENORK3Extend(obj, C, iteration)
		NewC = WENO5RK3Extend(obj, C, iteration)
		
		% extend C field away for surface of auxilary level set
		NewC = AENORK2Extend(obj, C, iter1, iter2, iter3)
		NewC = ANENORK2Extend(obj, C, iter1, iter2, iter3)

		% surface redistance schemes of various orders of accuracy
		NewA = ENORK2ClosetPointSurfaceRedistance(obj,A,iter1,iter2)
		NewA = WENORK3ClosetPointSurfaceRedistance(obj,A,iter1,iter2)
		NewA = WENO5RK3ClosetPointSurfaceRedistance(obj,A,iter1,iter2)
	end

	methods
		NewC = WENORK3SurfaceConservationLaw(obj,C,vx,vy,vz,iter,dt);
		NewC = SurfaceConservationLaw(obj,C,vx,vy,vz,iter,dt);
	end

	methods
		[fx,fy,fz] = GradientWENO(obj,f)
		[fxx,fyy,fzz,fxy,fyz,fzx] = HessianWENO(obj,f)
		val = LaplacianWENO(obj,f)
		[fx,fy,fz,fxx,fyy,fzz,fxy,fyz,fzx,f_lap] = GradientHessianLaplacian(obj,f)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization methods : 91 lines
	methods 

		% plot a 3D field on the val contour of the distance function
		function plotField(obj,val,Field,EdgeAlpha)
			% triangle mesh of the val isosurface. 
			% TriMesh is a structure with fields "vertices" and "faces"
			TriMesh = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,val);
			% interpolate the values of the Field onto the vertices of the triangle mesh
			SurfField = interp3(obj.GD3.X, obj.GD3.Y, obj.GD3.Z, Field, ...
				TriMesh.vertices(:,1), TriMesh.vertices(:,2), TriMesh.vertices(:,3), 'linear');
			% plot surface mesh 
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces, ...
				  'FaceVertexCData',SurfField,'FaceColor','interp',...
				  'EdgeColor','k','EdgeAlpha',EdgeAlpha)
			axis equal
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces,'FaceVertexCData',SurfField,...
				'FaceColor','interp','EdgeColor','none')
			axis equal
			view(45,30)
			colorbar
			camlight; lighting gouraud
		end

		% plot several half contours of the distance function
		function plot(obj)
			axis(obj.GD3.BOX)
			obj.plotIso(-12*obj.GD3.Dx,0.8,'Red')
			obj.plotIso(-6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(0,0.8,'Blue')
			obj.plotIso(6*obj.GD3.Dx,0.8,'Green')
			obj.plotIso(12*obj.GD3.Dx,0.8,'Red')
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end

		% plot half of the val contour of the distance function
		function plotIso(obj,val,trans,Color)
			F = obj.F;
			F(obj.GD3.Y<0) = inf;
			%F((obj.GD3.X+obj.GD3.Y)<0) = inf;
			%F( obj.GD3.X<0 & obj.GD3.Y<0) = inf;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
		end

		% plot the val contour of the distance function
		function plotSurface(obj,val,trans,FaceColor, EdgeColor)
			F = obj.F;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',FaceColor,'EdgeColor',EdgeColor,'FaceAlpha',trans,'EdgeAlpha',0.5);
			axis(obj.GD3.BOX)
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end

		% plot the val contour of the field (not the level set function)
		function plotSurfaceField(obj,F,val,trans,Color)
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
			axis(obj.GD3.BOX)
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end

		% plot the val contour of the field within the surface (designed for the auxilary level
	 	% set function
		[x,y,z]=CrossingLine(obj,iso,field, faces, verts, colors) % defined elsewhere
		function plotIsoField(obj, iso, field, PlotSurface)
			[faces,verts,colors] = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,0,field);
			if PlotSurface
				%obj.plotSurface(0,1,'Green','Black');
				obj.plotField(0,obj.AHeaviside,0.01)
			end
			[x,y,z] = obj.CrossingLine(0,field,faces,verts,colors);
			line(x(:),y(:),z(:),'Color','red','LineWidth',3);
			len = length(iso);
			for i=1:len
				[x,y,z] = obj.CrossingLine(iso(i),field,faces,verts,colors);
				line(x(:),y(:),z(:),'Color','black','LineWidth',2);
			end
		end

	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
