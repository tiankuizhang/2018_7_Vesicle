% test scheme for constructing cubic interpolants and root finding

%
s = 0.02;

v0 = 2; 
v1 = .8;
v2 = -0.1;
v3 = 1.2;

xl = [0, s, 2*s, 3*s];
yl = [v0, v1, v2, v3];


[c0, c1, c2, c3] = CubicInterpolant(v0,v1,v2,v3,s);

xc = CubicRoot(c0,c1,c2,c3,s,v1,v2)

fun = @(x) c0 + c1*x + c2*x.^2 + c3*x.^3;

x = 0:0.0001:3*s;
y = fun(x);

figure
plot(x,y,'Color','blue')
line(xl,yl,'Color','green')

line([0,3*s],[0,0],'Color','red')
line([xc,xc],[min(yl),max(yl)],'Color','red')






% given (v0,v1,v2,v3) such that v1*v2 <0, we want to find a cubic polynomial that
% interpolates through those four points
% denote it by q3  = c0 + c1*x + c2*x^2 + c3*x^3 ;

% let use x0 as the origin, then c0 = v0
% v1-v0 =   s*c_1 +   s^2*c_2 +    s^3*c_3   ( s  s^2   s^3 ) ( c_1 ) 
% v2-v0 = 2*s*c_1 + 4*s^2*c_2 +  8*s^3*c_3 = (2s 4s^2  8s^3 ) ( c_2 )
% v3-v0 = 3*s*c_1 + 9*s^2*c_2 + 27*s^3*c_3   (3s 9s^2 27s^3 ) ( c_3 )

%syms s
%A = [    s,   s^2,    s^3; ...
%	   2*s,	4*s^2,  8*s^3; ...
%	   3*s,	9*s^2, 27*s^3;] 
%
% inv(A)
% [		3/s,	-3/(2*s), 		1/(3*s)		]
% [-5/(2*s^2), 	2/s^2,			-1/(2*s^2)	]
% [1/(2*s^3),	-1/(2*s^3), 	1/(6*s^3)	]

% c_1 =      (3/s) * (v1-v0) - (3/(2*s)) * (v2-v0) + (1/(3*s)) * (v3-v0);
% c_2 = -5/(2*s^2) * (v1-v0) +     2/s^2 * (v2-v0) - 1/(2*s^2) * (v3-v0);
% c_3 =  1/(2*s^3) * (v1-v0) - 1/(2*s^3) * (v2-v0) + 1/(6*s^3) * (v3-v0);


% given four points (0,v0),(s,v1),(2*s,v2),(3*s,v3)
% return the coefficients of the cubic interpolant 
% denote it by q3  = c0 + c1*x + c2*x^2 + c3*x^3 ;
function [c0, c1, c2, c3] = CubicInterpolant(v0,v1,v2,v3,s)

	c0 = v0;
 	c1 =        3/s * (v1-v0) -   3/(2*s) * (v2-v0) +   1/(3*s) * (v3-v0);
 	c2 = -5/(2*s^2) * (v1-v0) +     2/s^2 * (v2-v0) - 1/(2*s^2) * (v3-v0);
 	c3 =  1/(2*s^3) * (v1-v0) - 1/(2*s^3) * (v2-v0) + 1/(6*s^3) * (v3-v0);
end

% find root between (s,v1) and (2s,v2)  with Newton's Method
function xc = CubicRoot(c0,c1,c2,c3,s,v1,v2)

	% initial guess from linear interpolant
	xc = s - s * v1 / (v2 - v1);

	iter = 0;
	max_iter = 50;
	epsilon = 1e-14;
	diff = 1;

	while diff>epsilon && iter < max_iter
		f = c0 + c1 * xc + c2 * xc^2 + c3 * xc^3;
		df = c1 + 2 * c2 * xc + 3 * c3 * xc^2;

		xc = xc - f / df;

		diff = abs(f / df);

		iter = iter + 1;
	end

	iter


end


