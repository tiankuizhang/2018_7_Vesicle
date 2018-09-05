classdef SDF3 < handle
	
	%SDF3 : signed distance function in 3D

	properties (SetAccess = immutable)
		GD3 % SD.GD3 object
	end

	properties
		F % values of the signed distance function
	end

	methods

		function obj = SDF3(grid, Xm, Ym, Zm, Val)
			obj.GD3 = grid;
			obj.GPUInitialize;
			obj.F = Val;
		end

	end

	%% GPU related properties and functions
	properties
		
		% parameter for GPU kernel
		ThreadBlockSize
		GridSize

		% kernel funcion object for ENORK2 reinitialization scheme
		ENORK2_boundary_correction % calculate grid step with modification near boundary
		ENORK2_reinitiliaztion_step % calculate the numerical Hamiltonian for the Reinitalization equation

		% kernel function object for ENORK2 extend scheme
		ENORK2_upwind_normal % calculate upwind normals of the level set function
		ENORK2_extend_step % calculate the extension step

		% kernel funcion object for ENORK2 surface redistance scheme
		ENORK2_surface_redistance_step % calculate the numerical Hamiltonian for the surface redistacne equation

	end

	methods
		
		% initialize GPU functions
		function GPUInitialize(obj)
			
			obj.ThreadBlockSize = [obj.GD3.mrows,ceil(512/obj.GD3.mrows),1];
			obj.GridSize = [ceil(obj.GD3.mrows/obj.ThreadBlockSize(1)), ...
							ceil(obj.GD3.ncols/obj.ThreadBlockSize(2)), ...
							ceil(obj.GD3.lshts/obj.ThreadBlockSize(3))];

			% functions used by reinitialization scheme and other schemes
			system('nvcc -ptx CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu -o CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx');

			obj.ENORK2_boundary_correction = parallel.gpu.CUDAKernel('CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx', ...
																	 'CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu', ...
																	 'boundary_correction');
			obj.ENORK2_boundary_correction.ThreadBlockSize = obj.ThreadBlockSize;
			obj.ENORK2_boundary_correction.GridSize = obj.GridSize;

			obj.ENORK2_reinitiliaztion_step = parallel.gpu.CUDAKernel('CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.ptx', ...
																	 'CUDA_Code/1_0_ENORK2_Reinitialization/boundary_correction.cu', ...
																	 're_step');
			obj.ENORK2_reinitiliaztion_step.ThreadBlockSize = obj.ThreadBlockSize;
			obj.ENORK2_reinitiliaztion_step.GridSize = obj.GridSize;

			% functions used by the extend scheme and other schemes
			system('nvcc -ptx CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu -o CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx');

			obj.ENORK2_upwind_normal = parallel.gpu.CUDAKernel('CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx', ...
				   											   'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu', ...
														   	   'upwind_normal');
	 		obj.ENORK2_upwind_normal.ThreadBlockSize = obj.ThreadBlockSize;
			obj.ENORK2_upwind_normal.GridSize = obj.GridSize;		

			obj.ENORK2_extend_step = parallel.gpu.CUDAKernel('CUDA_Code/2_0_ENORK2_Extend/enork2_extend.ptx', ...
				   											   'CUDA_Code/2_0_ENORK2_Extend/enork2_extend.cu', ...
														   	   'extend_step');
	 		obj.ENORK2_extend_step.ThreadBlockSize = obj.ThreadBlockSize;
			obj.ENORK2_extend_step.GridSize = obj.GridSize;		

			% functions used by the surface redistance schemes
			system('nvcc -ptx CUDA_Code/3_0_ENORK2_SurfaceRedistance/enork2_surface_redistance.cu -o CUDA_Code/3_0_ENORK2_SurfaceRedistance/enork2_surface_redistance.ptx');

			obj.ENORK2_surface_redistance_step = parallel.gpu.CUDAKernel('CUDA_Code/3_0_ENORK2_SurfaceRedistance/enork2_surface_redistance.ptx', ...
																		 'CUDA_Code/3_0_ENORK2_SurfaceRedistance/enork2_surface_redistance.cu', ...
																		 'surface_redistance_step');
			obj.ENORK2_surface_redistance_step.ThreadBlockSize = obj.ThreadBlockSize;
			obj.ENORK2_surface_redistance_step.GridSize = obj.GridSize;
		end
		

	end


	methods
		ENORK2Reinitialization(obj,iteration)	
		NewC = ENORK2Extend(obj, C, iteration)
		NewAF = ENORK2CentralUpwindSurfaceRedistance(obj,AF,iteration)
	end

	%%
	%% visualization methods
	%%
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
			F(obj.GD3.Y<0) = inf;
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
		function plotIsoField(obj, iso, field)
			[faces,verts,colors] = isosurface(obj.GD3.X,obj.GD3.Y,obj.GD3.Z,obj.F,0,field);
			obj.plotSurface(0,1,'Green',1);
			[x,y,z] = obj.CrossingLine(0,field,faces,verts,colors);
			line(x(:),y(:),z(:),'Color','red','LineWidth',3);
			len = length(iso);
			for i=1:len
				[x,y,z] = obj.CrossingLine(iso(i),field,faces,verts,colors);
				line(x(:),y(:),z(:),'Color','black','LineWidth',2);
			end
		end

		% calculate intersection of the level suface of obj.F and iso level curve of field
		function [x,y,z]=CrossingLine(obj,iso,field, faces, verts, colors)

			mask = colors>iso; % true if verts is outside

			outcount = sum(mask(faces),2); % 3 for outside, 0 for inside
			cross = (outcount == 2) | (outcount == 1);
			cross_tris = faces(cross,:);

			% make all cross_tris to be 1 vertex inside and 2 outside
			% since they can be treated in the same way
			out_vert = mask(cross_tris);
			flip = sum(out_vert,2) == 1;
			out_vert(flip,:) = 1-out_vert(flip,:);

			ntri = size(out_vert,1);
			overt = zeros(ntri,3);
		
			% now the first element is the only one outside/inside	
			for i=1:ntri
				v1i = find(~out_vert(i,:));
				v2i = 1 + mod(v1i,3);
				v3i = 1 + mod(v1i+1,3);
				overt(i,:) = cross_tris(i,[v1i v2i v3i]);
			end
			
			% value for linear interpolation	
			u = (iso - colors(overt(:,1))) ./ (colors(overt(:,2)) - colors(overt(:,1)));
			v = (iso - colors(overt(:,1))) ./ (colors(overt(:,3)) - colors(overt(:,1)));
			
			% convert linear interpolation values to x,y,z coordinates
			uverts = repmat((1-u),[1 3]).*verts(overt(:,1),:) + ...
					 repmat(    u,[1 3]).*verts(overt(:,2),:);
			vverts = repmat((1-v),[1 3]).*verts(overt(:,1),:) + ...
				     repmat(v    ,[1 3]).*verts(overt(:,3),:);
			
			% construct line segments
			% to plot: line(x(:),y(:),z(:),'Color','Red','LineWidth',3)
			x = nan(3,ntri);
			x(1,:) = uverts(:,1)';
			x(2,:) = vverts(:,1)';
			y = nan(3,ntri);
			y(1,:) = uverts(:,2)';
			y(2,:) = vverts(:,2)';
			z = nan(3,ntri);
			z(1,:) = uverts(:,3)';
			z(2,:) = vverts(:,3)';

		end


	end

end
