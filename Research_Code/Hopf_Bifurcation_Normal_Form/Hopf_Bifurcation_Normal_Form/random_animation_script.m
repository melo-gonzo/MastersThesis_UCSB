% animation script 
% run main first to get xgamma and xp
close all
vid=1;
step_size=20;
pause_length=0.001;
xptheta=a_phase(p,sqrt(s.xp_short(end,1)^2+s.xp_short(end,2)^2),atan2(s.xp_short(end,2),s.xp_short(end,1)));
xtheta=a_phase(p,sqrt(s.x_short(end,1)^2+s.x_short(end,2)^2),atan2(s.x_short(end,2),s.x_short(end,1)));
rr=0.001:0.001:1;
r2=sqrt(s.xp_short(end,2)^2+s.xp_short(end,1)^2);
xpiso=xptheta+(p.d(1)/p.c(1))*log(rr)-(p.d(1)/(2*p.c(1)))*log(-p.a(1)/p.c(1));
xiso=xtheta+(p.d(1)/p.c(1))*log(rr)-(p.d(1)/(2*p.c(1)))*log(-p.a(1)/p.c(1));

%% Plot Phase Field
figure(); hold on;
[xm,ym]=meshgrid(-1:0.001:1,-1:0.001:1);
[m,n]=size(xm);
theta=nan(m,n);
for ii=1:m
    for jj=1:n
        theta(ii,jj)=a_phase(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2),...
            atan2(ym(ii,jj),xm(ii,jj)));
        %           psi_field(ii,jj)=a_approach(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2));
    end
end
pcolor(xm,ym,mod(theta,2*pi));shading flat
colormap(flip([cool;flipud(cool)]));axis square;box on
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
clear theta xm ym n m ii jj


%% plot xgamma
if vid==1
    v=VideoWriter('hopf_simulation_iso_straight');
    open(v)
end
plot(s.xg(:,1),s.xg(:,2),'k','linewidth',1);
title('Periodic Orbit')
if vid==1;writeVideo(v,getframe(figure(1)));end
%% Animate one rotation of the periodic orbit

xgp=plot(s.xg(1,1),s.xg(1,2),'b.','markersize',30);
for k=1:step_size:length(s.xg(:,1))
    set(xgp,'XData',s.xg(k,1),'YData',s.xg(k,2));
    if vid==1;writeVideo(v,getframe(figure(1)));end
    pause(pause_length)
end

del_last_line(1)
%% Plot periodic and perturbed
xgp = plot(s.x(1,1),s.x(1,2),'b.','markersize',30);
xpp = plot(s.xp(1,1),s.xp(1,2),'r.','markersize',30);
title('Periodic and Perturbed Orbit')
if vid==1;writeVideo(v,getframe(figure(1)));end
for k=1:step_size:2*length(s.xg(:,1))
    set(xgp,'XData',s.x(k,1),'YData',s.x(k,2));
    set(xpp,'XData',s.xp(k,1),'YData',s.xp(k,2));
    if k==p.tau_idx
        line_flash=[0.1:0.005:3 3:-0.005:0.1 0.1:0.005:3 3:-0.005:0.1];
        iso1=plot(rr.*cos(xpiso),rr.*sin(xpiso),'r-.','linewidth',0.1);
        iso2=plot(rr.*cos(xiso),rr.*sin(xiso),'b-.','linewidth',0.1);
        title('Phase Change')
        for n=1:step_size:length(line_flash)
            set(iso1,'linewidth',line_flash(n));
            set(iso2,'linewidth',line_flash(n));
            if vid==1;writeVideo(v,getframe(figure(1)));end
            pause(pause_length)
        end
        del_last_line(2)
        psi1=plot(p.x0(1)*cos(0:0.001:2*pi),p.x0(1)*sin(0:0.001:2*pi),'b-.',...
            'linewidth',0.1);
        psi2=plot(r2*cos(0:0.001:2*pi),r2*sin(0:0.001:2*pi),'r-.',...
            'linewidth',0.1);
        title('Transient Amplitude Change')
        for n=1:step_size:length(line_flash)
            set(psi1,'linewidth',line_flash(n));
            set(psi2,'linewidth',line_flash(n));
            if vid==1;writeVideo(v,getframe(figure(1)));end
            pause(pause_length)
        end
        del_last_line(2)
        title('Periodic and Perturbed Orbit')
    end
    if vid==1;writeVideo(v,getframe(figure(1)));end
    pause(pause_length)
end
if vid==1
    close(v)
end







