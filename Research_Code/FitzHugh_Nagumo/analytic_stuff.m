% close all
% clear all
% clc

c=-1;
d=1;
b=1;
a=0.05;
a2=1.1;
c1=1000;
t=0:0.001:2*pi;
theta=0:0.001:2*pi;

r1=1i*sqrt(a)*exp(a*(c1+t))./sqrt(-1+c*exp(2*a*(c1+t)));
c2=-2*pi;
r2=1i*sqrt(a2)*exp(a2*(c2+t))./sqrt(-1+c*exp(2*a2*(c2+t)));
% % 
r=r2+r1(end);
% 
x=r.*cos(theta);
y=r.*sin(theta);

plot(x,y,'r','linewidth',2);


%%
% p=parameters(0.05,1,-1,1);
% p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
% p.dt=0.001;
% p.tspan=0:p.dt:30;
% p.tau=1;
% p.dp=p.a+0.1;
% p.a=[p.a p.a+p.dp p.a];
% 
% [t,tp,x,xp]=hopf_ode(p);
% 
% hold on;plot(xp(1:1000,1),xp(1:1000,2),'r')
% plot(x(:,1),x(:,2),'k')