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


	%%
	%% GPU related properties and functions
	%%
	properties
		
		% parameter for GPU kernel
		ThreadBlockSize
		GridSize

		% kernel funcion object for ENORK2 reinitialization scheme
		ENORK2_boundary_correction % calculate grid step with modification near boundary
		ENORK2_reinitiliaztion_step % calculate the numerical Hamiltonian for the Reinitalization equation

	end

	methods
		
		% initialize GPU functions
		function GPUInitialize(obj)
			
			obj.ThreadBlockSize = [obj.GD3.mrows,4,1];
			obj.GridSize = [ceil(obj.GD3.mrows/obj.ThreadBlockSize(1)), ...
							ceil(obj.GD3.ncols/obj.ThreadBlockSize(2)), ...
							ceil(obj.GD3.lshts/obj.ThreadBlockSize(3))];

			% functions used by reinitialization scheme
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

		end
		

	end


	methods
		
		function ENORK2Reinitialization(obj)
			
			Fgpu = obj.F;
			
			mask = obj.F<0;
			deltat = zeros(obj.GD3.Size, 'gpuArray');

			xpr = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
			xpl = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dx;
			ypf = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
			ypb = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dy;
			zpu = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;
			zpd = ones(obj.GD3.Size, 'gpuArray') * obj.GD3.Dz;

			[xpr, xpl, ypf, ypb, zpu, zpd] = feval(obj.ENORK2_boundary_correction, ...
					xpr, xpl, ypf, ypb, zpu, zpd, Fgpu, obj.GD3.NumElt, ...
				    obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
					obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz);	

			minx = min(xpr,xpl);
			miny = min(ypf,ypb);
			minz = min(zpu,zpd);
			deltat = 0.3 * min(minx, min(miny,minz));

			step = zeros(obj.GD3.Size, 'gpuArray');

			for i=1:100
				step = feval(obj.ENORK2_reinitiliaztion_step, step, Fgpu, mask, deltat, ...
						xpr, xpl, ypf, ypb, zpu, zpd, ...
				    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
						obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
				Ftmp = Fgpu - step;
				step = feval(obj.ENORK2_reinitiliaztion_step, step, Ftmp, mask, deltat, ...
						xpr, xpl, ypf, ypb, zpu, zpd, ...
				    	obj.GD3.mrows, obj.GD3.ncols, obj.GD3.lshts, ...
						obj.GD3.Dx, obj.GD3.Dy, obj.GD3.Dz, obj.GD3.NumElt);	
				Fgpu = (Fgpu + Ftmp - step) / 2;
			end

			obj.F = Fgpu;

		end

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
			set(p1,'FaceColor',Color,'EdgeColor','k','FaceAlpha',trans);
			axis(obj.GD3.BOX)
			daspect([1 1 1])
			view(3); 
			camlight; lighting gouraud
		end
	end

end
