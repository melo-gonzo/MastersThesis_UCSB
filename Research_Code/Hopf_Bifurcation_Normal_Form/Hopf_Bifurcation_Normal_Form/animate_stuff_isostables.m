close all
clear all
beep off
clc
set_default_plot

p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*4,3);
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
        theta_field(ii,jj)=a_approach(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2));
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
r=0.001:0.001:1.55*max(s.xg(:,1));
r=r(end)*ones(length(r),1);
theta=-a_approach(p,r);
% plot(r.*cos(theta),r.*sin(theta),'w-.','linewidth',2)

%% Plot points on isochron
n_pts=10;
x0_idx=[(round(1:length(theta)/(n_pts):length(theta))) length(theta)];
x0_idx(1)=[];
tts=0:(2*pi)/length(x0_idx):2*(pi);
tts(end)=[];
x0=[(r(x0_idx)).*cos(tts)' (r(x0_idx)).*sin(tts)'];
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
    v=VideoWriter('isostable_simulation');
    open(v)
end


pcolor(xm,ym,exp(theta_field));shading flat
colormap(jet)
axis square; axis tight;xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.125)
% clear xm ym n m ii jj
set(gca,'YTick',[-0.5 0 0.5]);
cb=colorbar('TickLabelInterpreter', 'latex');
set(get(cb,'title'),'string','$e^{\psi}$','interpreter','latex');
line_flash=[0.1:0.005:3 3:-0.005:0.1 0.1:0.005:3 3:-0.005:1.5];
if vid==1;writeVideo(v,getframe(figure(1)));end

step_size=20;
hold on;xg_plot=plot(s.xg(:,1),s.xg(:,2),'w'); title('Periodic Trajectory')
for n=1:step_size:length(line_flash)
    set(xg_plot,'linewidth',line_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 
set(xg_plot,'linewidth',1.5)
title('Isostable Level Set')
isochron=plot(r(1).*cos(0:0.01:2*pi),r(1).*sin(0:0.01:2*pi),'w-.','linewidth',2)
for n=1:step_size:length(line_flash)
    set(isochron,'linewidth',line_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 
set(isochron,'linewidth',1.5)

title('Initial Conditions')
lines=plot(x0(:,1),x0(:,2),'k.','markersize',20);
axis tight
m_flash=line_flash*10;
for n=1:step_size:length(m_flash)
    set(lines,'markersize',m_flash(n))
    if vid==1;writeVideo(v,getframe(figure(1)));end
end 

set(lines,'markersize',20)
title('Flow')
for k=1:step_size:length(p.tspan)
    set(lines,'XData',s_x(k,:),'YData',s_y(k,:))
    if vid==1;writeVideo(v,getframe(figure(1)));end
    pause(0.01)
end 
if vid==1
    close(v)
end





