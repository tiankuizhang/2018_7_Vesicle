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

