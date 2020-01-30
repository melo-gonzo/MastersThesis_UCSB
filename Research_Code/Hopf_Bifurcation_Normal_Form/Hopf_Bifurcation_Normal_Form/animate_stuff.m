close all
clear all
beep off
clc
set_default_plot

p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*8,3);
numeric=1;
p.perturbed=0;  
pert=[1 0];
[s.xg, s.tg]=xg_function(p,numeric);

%% Plot Phase Field 
[xm,ym]=meshgrid(-0.5:0.001:0.5,-0.5:0.001:0.5);
[m,n]=size(xm);
theta_field=nan(m,n);
for ii=1:m
    for jj=1:n
        theta_field(ii,jj)=a_phase(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2),...
            atan2(ym(ii,jj),xm(ii,jj)));
    end
end 
% pcolor(xm,ym,mod(theta_field,2*pi));shading flat
% colormap(jet);
% axis square; axis tight;xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.125)
% clear theta xm ym n m ii jj
% set(gca,'YTick',[-0.5 0 0.5])
%% Plot xg
% hold on;plot(s.xg(:,1),s.xg(:,2),'w')

%% Plot Single Isochron
r=0.001:0.001:1.75*max(s.xg(:,1));
theta=-a_phase(p,r,-5*pi/12);
% plot(r.*cos(theta),r.*sin(theta),'w-.','linewidth',2)

%% Plot points on isochron
n_pts=10;
x0_idx=[(round(1:length(theta)/(n_pts):length(theta))) length(theta)];
x0_idx(1)=[];
x0=[(r(x0_idx).*cos(theta(x0_idx)))' (r(x0_idx).*sin(theta(x0_idx)))'];
% lines=plot(x0(:,1),x0(:,2),'k.','markersize',20);
s_x=zeros(length(p.tspan),length(x0_idx));
s_y=zeros(length(p.tspan),length(x0_idx));
for k=1:length(x0_idx)
    p.x0=x0(k,:);
    [~,s.tp,s.x,s.xp]=hopf_ode(p);
    s_x(:,k)=s.x(:,1);
    s_y(:,k)=s.x(:,2);
end 


%% Make Animation
vid=1;
if vid==1
    v=VideoWriter('isochron_simulation');
    open(v)
end


pcolor(xm,ym,mod(theta_field,2*pi));shading flat
colormap(jet);
axis square; axis tight;xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.125)
% clear xm ym n m ii jj
set(gca,'YTick',[-0.5 0 0.5])
cb=colorbar('TickLabelInterpreter', 'latex');set(cb,'Ticks',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
set(cb,'TickLabels',xtl)
set(get(cb,'title'),'string','$\theta$','interpreter','latex');
line_flash=[0.1:0.005:3 3:-0.005:0.1 0.1:0.005:3 3:-0.005:1.5];
if vid==1;writeVideo(v,getframe(figure(1)));end

step_size=20;
hold on;xg_plot=plot(s.xg(:,1),s.xg(:,2),'w'); title('Periodic Trajectory')
for n=1:step_size:length(line_flash)
    set(xg_plot,'linewidth',line_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 
set(xg_plot,'linewidth',1.5)
title('Isochron Level Set')
isochron=plot(r.*cos(theta),r.*sin(theta),'w-.','linewidth',2)
for n=1:step_size:length(line_flash)
    set(isochron,'linewidth',line_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 
set(isochron,'linewidth',1.5)

title('Initial Conditions')
lines=plot(x0(:,1),x0(:,2),'k.','markersize',20);
m_flash=line_flash*10;
for n=1:step_size:length(m_flash)
    set(lines,'markersize',m_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 

set(lines,'markersize',20)
title('Flow')
for k=1:step_size*2:length(p.tspan)
    set(lines,'XData',s_x(k,:),'YData',s_y(k,:))
    if vid==1;writeVideo(v,getframe(figure(1)));end
    pause(0.01)
end 
if vid==1
    close(v)
end





