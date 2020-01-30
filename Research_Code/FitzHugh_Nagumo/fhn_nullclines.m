function fhn_nullclines(a,b,c1,c2,d)
x=-2.5:0.001:2.5;
y=-2.5:0.001:2.5;


xn=x-x.^3;
yn=(x+a)/b;

plot(x,xn);hold on;plot(y,yn);axis tight 

% figure(1);
[xp,yp]=meshgrid(-2.5:0.2:2.5,-2.5:0.2:2.5);
xv=xp-xp.^3 -yp;
yv=d*(xp+a-b*yp);
quiver(xp,yp,xv,yv,'k-','linewidth',0.5);axis tight
axis([-2.5 2.5 -2.5 2.5])