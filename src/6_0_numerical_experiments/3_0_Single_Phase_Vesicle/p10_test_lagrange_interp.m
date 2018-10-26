
N = 4;
x = linspace(0,1,N)
y = sin(x)

kk = zeros(1,N);
xx = zeros(N,N-1);

for i=1:N
	kk(i) = y(i);
	ind = 1;
	for j=1:N
		if j~=i
			kk(i) = kk(i) / (x(i) - x(j));
			xx(i,ind) = x(j);
			ind = ind + 1;
		end
	end
end

c = zeros(N,N);
c(:,1) = 1;

for i=1:N
	for j=1:N-1
		c(i,2:(j+1)) = c(i,2:(j+1)) - xx(i,j) .* c(i,1:j);
	end
end

coef = zeros(1,N);

for i=1:N
	for j=1:N
		coef(i) = coef(i) + kk(j) * c(j,N+1-i);
	end
end
