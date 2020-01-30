close all
clear all
beep off
clc
set_default_plot
%% Parameters Setup
close all
x=-2.5:0.001:2.5;
y=-2.5:0.001:2.5;
xn=x-x.^3;
line2=plot(0,0,'b-');hold on;line3=plot(0,0,'b-');
[xp,yp]=meshgrid(-2.5:0.2:2.5,-2.5:0.2:2.5);


epss=0.8;
line=plot(0,0);
axis([-2 2 -2 2])

for a=0.09:0.01:0.09
    p=parameters(-a,0.8,0,0,epss);
    p.x0=[-1 0];
    p.dt=0.01;
    p.tspan=0:p.dt:5000;
    p.perturbed=0;
    
    yn=(x+p.a)/p.b;
    yv=p.d*(xp+p.a-p.b*yp);
    set(line2,'XData',x,'YData',xn);
    set(line3,'XData',y,'YData',yn);
    
    [s.t,~,s.x,~]=fhn_ode(p);
    set(line,'XData',s.x(:,1),'YData',s.x(:,2));
    pause(0.001)
    title(num2str(a))
end







% fhn_nullclines(0.1,0.8,0,0,epss);
% plot(s.x(:,1),s.x(:,2))