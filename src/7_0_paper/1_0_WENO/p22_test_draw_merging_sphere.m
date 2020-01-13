load('mergingSphere')

figure

subplot(1,2,1)
r = 0.6;
semilogy(Dist./r,EM(1,:),'-*')
hold on
semilogy(Dist./r,EM(2,:),'-x')
semilogy(Dist./r,EM(3,:),'-o')
semilogy(Dist./r,EM(4,:),'-s')
grid on
legend( ' 16\times 16\times 32', ...
		' 32\times 32\times 64', ...
		' 64\times 64\times128', ...
		'128\times128\times256', ...
		'Location', 'southwest')
xlabel('distacne between two spheres')
ylabel('||\phi - \phi_h||_{\infty}')
title('deterioration of accuracy in the presence of kinks')

subplot(1,2,2)
dx = 2 ./ [15,31,63,127];
index = [19,3];
loglog(dx, EM(:,19),'-*')
hold on
loglog(dx, EM(:,13),'-x')
loglog(dx, EM(:,12),'-o')
loglog(dx, EM(:,10),'-s')
grid on
legend( [sprintf('%.1f', Dist(19)/r) + "\timesr"], ...
		[sprintf('%.1f', Dist(13)/r) + "\timesr"], ...
		[sprintf('%.1f', Dist(12)/r) + "\timesr"], ...
		[sprintf('%.1f', Dist(10)/r) + "\timesr"], ...
		'Location', 'southeast')
xlabel('grid spacing')
ylabel('||\phi - \phi_h||_{\infty}')
title('deteroriation of convergence rate in the presence of kinks ')

