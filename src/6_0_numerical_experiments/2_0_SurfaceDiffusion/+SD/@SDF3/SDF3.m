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
	end


	methods
		% update the above properties
		setCalculustToolBox(obj)
		GPUsetCalculusToolBox(obj)  
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
		AFDiracDelta % product of ADiracDelta and FDiracDelta

	end

	methods
		AsetCalculusToolBox(obj)
		GPUAsetCalculusToolBox(obj)
	end

	methods
		% return length of curve
		function length = calLength(obj)
			mask = abs(obj.F) < 2*obj.GD3.Ds & abs(obj.A) < 2*obj.GD3.Ds;
			length = sum(obj.AFDiracDelta(mask).*obj.NormCrossAF(mask)) * obj.GD3.Ds.^3;
		end
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GPU related properties and functions: 27 lines
	properties
		
		% parameter for GPU kernel
			ThreadBlockSize
			GridSize

		% kernel funcions object for ENORK2 reinitialization scheme
			% calculate grid step with modification near boundary
			ENORK2_boundary_correction 
			% calculate the numerical Hamiltonian for the Reinitalization equation
			ENORK2_reinitiliaztion_step 

		% kern functio oect forer WENORK3 reinitialization  scheme
			% calculate distance to the nearby nodes including boundary nodes
			cubic_boundary_correction
			% calculate the reinitialization step
			WENORK3_reinitialization_step


		% kernel functions object for ENORK2 extend scheme
			ENORK2_upwind_normal % calculate upwind normals of the level set function
			ENORK2_boundary_interpolate % interpolate values at the boundary
			ENORK2_extend_step % calculate the extension step

		% kernel functions object for WENORK3 extend scheme
			WENORK3_upwind_normal % calculate upwind normals of the level set function
			WENORK3_boundary_interpolate % interpolate values at the boundary
			WENORK3_extend_step % calculate the extension step

		% kernel function object for GPUsetCalculusToolBox scheme
			set_calculus_toolbox % set Fx,Fy ...
			auxi_set_calculus_toolbox % set Ax,Ay ...

	end
	
	methods 
		GPUInitialize(obj)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% utilities : reinitliazation, extend, surface_redistance
	methods
		NewF = ENORK2Reinitialization(obj,F,iteration)	
		NewF = WENORK3Reinitialization(obj,F,iteration)	
		NewC = ENORK2Extend(obj, C, iteration)
		NewC = WENORK3Extend(obj, C, iteration)
		NewA = ENORK2ClosetPointSurfaceRedistance(obj,A,iter1,iter2)
		NewA = WENORK3ClosetPointSurfaceRedistance(obj,A,iter1,iter2)
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualization methods : 86 lines
	methods 

		% plot a 3D field on the val contour of the distance function
		function plotField(obj,val,Field)
			% triangle mesh of the val isosurface. 
			% TriMesh is a structure with fields "vertices" and "faces"
			TriMesh = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,val);
			% interpolate the values of the Field onto the vertices of the triangle mesh
			SurfField = interp3(obj.GD3.X, obj.GD3.Y, obj.GD3.Z, Field, ...
				TriMesh.vertices(:,1), TriMesh.vertices(:,2), TriMesh.vertices(:,3), 'linear');
			% plot surface mesh 
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces, ...
				  'FaceVertexCData',SurfField,'FaceColor','interp','EdgeColor','none')
			axis equal
			patch('Vertices',TriMesh.vertices,'Faces',TriMesh.faces,'FaceVertexCData',SurfField,...
				'FaceColor','interp','EdgeColor','none')
			axis equal
			view(45,30)
			colorbar
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
			%F(obj.GD3.Y<0) = inf;
			F((obj.GD3.X+obj.GD3.Y)<0) = inf;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
		end

		% plot the val contour of the distance function
		function plotSurface(obj,val,trans,Color, time)
			F = obj.F;
			surf1 = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,val);
			p1 = patch(surf1);
			isonormals(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,F,p1)
			set(p1,'FaceColor',Color,'EdgeColor','none','FaceAlpha',trans);
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
				obj.plotSurface(0,1,'Green',1);
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
