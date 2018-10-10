% MeshFree Approximation Methods with Matlab Program1.1.
function DM = DistanceMatrix(dsites,ctrs)

	[M,s] = size(dsites); [N,s] = size(ctrs);

	DM = zeros(M,N);

	for d=1:s
		[dr,cc] = ndgrid(dsites(:,d),ctrs(:,d));
		DM = DM + (dr-cc).^2;
	end

	DM = sqrt(DM);
		
