close all
clear all
beep off
clc
set_default_plot
%% Parameters Setup
line1=plot3(0,0,0);
for ib=0:1:50
    p=parameters(1,0.05,-70,3,50,5,-90,5,0,ib);
    TPERT=0;
    p.x0=[-6.6507 0.24773 0.0017566];
    p.dt=0.01;
    p.tspan=0:p.dt:500;
    p.perturbed=0;
    [s.t,~,s.x,~]=thalamic_ode(p);
    set(line1,'XData',s.x(:,1),'YData',s.x(:,2),'ZData',s.x(:,3))
    title(strcat('ib: ',num2str(ib)))
    pause(0.001)
end







% fhn_nullclines(0.1,0.8,0,0,epss);
% plot(s.x(:,1),s.x(:,2))