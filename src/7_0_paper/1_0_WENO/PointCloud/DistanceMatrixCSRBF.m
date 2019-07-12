function DM = DistanceMatrixCSRBF(dsites,ctrs,ep)

	N = size(dsites,1); M = size(ctrs,1);

	support = 1.0/ep;
	nzmax = 25*N; 
	rowidx = zeros(1,nzmax); colidx = zeros(1,nzmax); validx = zeros(1,nzmax);
	istart = 1; iend = 0;
	if M > N
		[tmp,tmp,Tree] = kdtree(ctrs,[]);
		for i = 1:N
			[pts,dist,idx] = kdrangequery(Tree,dsites(i,:),support);
			newentries = length(idx);
			iend = iend + newentries;
			rowidx(istart:iend) = repmat(i,1,newentries);
			colidx(istart:iend) = idx';
			validx(istart:iend) = 1-ep*dist';
			istart = istart + newentries;
		end
	else
		[tmp,tmp,Tree] = kdtree(dsites,[]);
		for j = 1:M
			[pts,dist,idx] = kdrangequery(Tree,ctrs(j,:),support);
			newentries = length(idx);
			iend = iend + newentries;
			rowidx(istart:iend) = idx';
			colidx(istart:iend) = repmat(j,1,newentries);
			validx(istart:iend) = 1-ep*dist';
			istart = istart + newentries;
		end
	end
	idx = find(rowidx);
	DM = sparse(rowidx(idx),colidx(idx),validx(idx),N,M);
	kdtree([],[],Tree);


