% locate boundary position with high order interpolant

Order = 5;

stencil = 1:Order;

% create coefficient matix
syms s
C = repmat( repmat(s, [Order,1]) .* stencil', [1,Order]) ...
	.^ repmat(stencil, [Order,1]);
invC = inv(C);
invC = matlabFunction(inv(C));

% 
fun = @(x) exp(x)-exp(pi);
N = 16;
x = linspace(2,4,N);
dx = x(2) - x(1);
y = fun(x);
boundaryNode = [];
for ii=1:N-1
	if y(ii)*y(ii+1)<0
		boundaryNode = [boundaryNode,ii];
	end
end

idx = boundaryNode(1);
stencil_idx = [idx-2,idx-1,idx,idx+1,idx+2,idx+3];
stencil_idx = min(N,max(1,stencil_idx))
stencil_x = x(stencil_idx) - x(stencil_idx(1));
stencil_y = y(stencil_idx);

% calculate the coefficients for the interpolation polynomial
c0 = stencil_y(1);
co = [c0; invC(dx) * (stencil_y(2:end)-stencil_y(1))'];
p = @(r) polyval(co(end:-1:1),r-x(stencil_idx(1)));
x0 = [x(idx),x(idx+1)];
%options = optimset('Display','iter');
%[numericalroot fval exitflag output] = fzero(p,x0,options)
numericalroot = zeroin(p,x0(1),x0(2),0)
relative_error = abs(numericalroot - pi) / pi

X = linspace(x(1),x(end),100);

subplot(1,2,1)
hold on
plot(X,fun(X))
plot(X,p(X))
scatter(x(stencil_idx),y(stencil_idx),50)
grid on
hold off

subplot(1,2,2)
hold on
plot(X,p(X)-fun(X))
scatter(x(stencil_idx),0*y(stencil_idx),50)
grid on
hold off

% a less complicated implementation of fzero based on Forsythe's 
% 'computer methods for mathematical computations'
% Using Brent's method, return the root of a function or functor know to lie between x1 and x2
% the root will be refined until its accuracy is tol
% from Numerical Recipes section 9.3
function b = zeroin(func,x1,x2,tol)

	ITMAX = 100; % maximum allowed number of iterations
	EPS = eps;

	% initialization
	a=x1; b=x2; c=x2;
	fa=func(a); fb=func(b); fc = fb;

	for ii=0:ITMAX
		% loop invariant: result should be bound by b,c; b is the best result; a is the previous result
		% if b,c is on the same side, exchange a,c so that result is bound by b,c
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0))
			c = a; fc = fa;
			d = b - a; e = d;
		end
		% make sure b is the best result
		if abs(fc) < abs(fb)
			a = b; b= c; c = a;
			fa = fb; fb = fc; fc = fa;
		end

		% check convergence
		% toler seems a better choice than the matlab default one
		toler = 2.0 * EPS * abs(b) + 0.5 * tol;
		%toler = 2.0 * tol * max(abs(b),1.0);
		xm = 0.5*(c-b);
		if (abs(xm) <= toler) || (fb == 0.0)
			break;
		end

		% choose bisection or interpolation
		if (abs(e) < toler) || (abs(fa) <= abs(fb))
			% bounds decreasing too slowly, use bisection
			d = xm; e = xm;
		else
			% attempt interpolation
			s = fb/fa;
			if (a == c)
				% linear interpolation 
				p = 2.0 * xm * s;
				q = 1.0 - s;
			else
				% inverse quadratic interpolation
				q = fa/fc;
				r = fb/fc;
				p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0));
				q = (q - 1.0)*(r - 1.0)*(s - 1.0);
			end
			if p > 0 % check wether in bounds
				q = -q;
			else
				p = - p;
			end
			% is interpolation acceptable
			min1 = 3.0 * xm * q - abs(toler*q);
			min2 = abs(e*q);
			if 2.0*p < min(min1,min2)
				% accept interpolation
				e = d; d = p/q;
			else
				% use bisection
				d = xm; e = d;
			end
		end

		% move last best guest to a
		a = b; 
		fa = fb;
		if abs(d) > toler % evaluate new trial root
			b = b + d;
		elseif b > c
			b = b - toler;
		else
			b = b + toler;
		end
		fb = func(b);
	end

	if ii == ITMAX
		display('Maximum number of iteration exceeded zbrent')
	end

	fprintf('%d iteraion used, feval: %.3e\n',ii,fb)
end
























