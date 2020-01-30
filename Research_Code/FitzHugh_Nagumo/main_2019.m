close all
clear all
beep off
clc
set_default_plot
%% Parameters Setup
% close all
epss=0.50;
ptrb=[0 0.01];
p=parameters(0,0.8,0,0,epss);
% TPERT=0;

p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);

p2=parameters(p.a+p.a*ptrb(1)+ptrb(2),0.8,0,0,epss);...
    p2.x0=[-1 0];p2.dt=0.001;p2.tspan=0:p2.dt:100;p2.perturbed=0;
[s.t2,~,s.x2,~]=fhn_ode(p2);
[tg1,xg1]=periodic_trajectory(s.t,s.x);
[tg2,xg2]=periodic_trajectory(s.t2,s.x2);
% fhn_nullclines(p.a+p.a*ptrb(1)+ptrb(2),0.8,0,0,epss);
% hold on
% plot(xg1(:,1),xg1(:,2),'g-')
% plot(xg2(:,1),xg2(:,2),'r-')
% axis([-1.5 1.5 -1.5 1.5])
% pause(0.001)
%% 
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;

%% Do 2d stuff
% plots2d
% figure();
% er=pcolor(a_vals,tau_vals,abs(final_vals));colormap jet;
% set(er,'EdgeColor','none')

%% Perturbation 
p.x0=s.xg(1,:);
p.tau=round(s.tg(end)*0.2,3);
p.tpert=round(s.tg(end)*3*pi/2/(2*pi),3);
% p.tpert=TPERT;
p.tpert_idx=round(p.tpert/p.dt,5)+1;
p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;

if p.tau_idx>(s.tg(end)/p.dt)
    s.prc=[s.prc;[s.prc(:,1)+s.prc(end,1) s.prc(:,2:end)]];
    s.irc=[s.irc;[s.irc(:,1)+s.irc(end,1) s.irc(:,2:end)]];
end 

p.perturbed=1;
p.dp=p.a(1)*ptrb(1)+ptrb(2);
p.p_time=p.tpert:p.dt:p.tpert+p.tau;
p.a=[p.a p.a+p.dp p.a];
p.phi_span=2*pi*(p.tpert:p.dt:p.tpert+p.tau)/s.tg(end);
[s.t_short,s.tp_short,s.x_short,s.xp_short]=fhn_ode(p);


s.t_short=s.t_short(p.tpert_idx:end);
s.tp_short=s.tp_short(p.tpert_idx:end);
s.x_short=s.x_short(p.tpert_idx:end,:);
s.xp_short=s.xp_short(p.tpert_idx:end,:);


p.perturbed=0;
[~,s.tp,s.x,s.xp]=fhn_ode(p);

%% Sensitivity Function
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end

%% Coordinate Recovery 
s.xy_cr=cr_function(p,s,0);
s.xy_recovered = s.x_short + s.xy_cr;
% 
% figure();
% subplot(3,1,1)
% hold on
% title('Trajectorys')
% plot(xg1(:,1),xg1(:,2),'k-')
% plot(xg2(:,1),xg2(:,2),'k-.')
% plot(s.x_short(:,1),s.x_short(:,2),'b-')
% plot(s.xp_short(:,1),s.xp_short(:,2),'r-')
% plot(s.xy_recovered(:,1),s.xy_recovered(:,2),'g-')
% axis tight
% subplot(3,1,2)
% hold on
% title('PRC')
% plot(s.tg,s.prc_og(:,2:end))
% plot(s.tg(p.tpert_idx),s.prc_og(p.tpert_idx,2:end),'r.','markersize',20)
% subplot(3,1,3)
% hold on
% title('IRC')
% plot(s.tg,s.irc_og(:,2:end))
% plot(s.tg(p.tpert_idx),s.irc_og(p.tpert_idx,2:end),'r.','markersize',20)
% get(gcf,'Position')
% set(gcf,'Position',[ans(1)+ans(3)/2 2 ans(3)/2 922.5])



%%%%%%%
%% Plot Coordinate Recovery Plots
figure(3);hold on
plot(s.xg(:,1),s.xg(:,2))
plot(s.x_short(:,1),s.x_short(:,2))
plot(s.xp_short(:,1),s.xp_short(:,2))
%% Coordinate Recovery Plot
plot(s.xy_recovered(:,1),s.xy_recovered(:,2))
%%%%%%%
%%
% legend('x_g','x_s','x_p','x_r')
% box on; 
% mval=max(max(abs(s.xy_recovered)));
% axis(1.3*[-mval mval -mval mval])
% axis square




% plot_stuff
%% Plot Phase Field 
% [xm,ym]=meshgrid(-1:0.001:1,-1:0.001:1);
% [m,n]=size(xm);
% theta=nan(m,n);
% for ii=1:m
%     for jj=1:n
%         if sqrt(xm(ii,jj)^2+ym(ii,jj)^2)>0.1
%         theta(jj,ii)=a_phase(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2),...
%             atan2(ym(ii,jj),xm(ii,jj)));
%         end
%     end
% end 
% pcolor(xm,ym,mod(theta,2*pi));shading flat
% colormap([jet;flipud(jet)]) 
% clear theta xm ym n m ii jj

%%
% figure()
% load final_coords_tpert_change_FHN_020805
% d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
% a_vals=0:0.003:1;
% tau_vals=round((0:0.003:1)*s.tg(end),3);
% dd=pcolor(a_vals,tau_vals,d)
% set(dd,'EdgeColor','none');colormap jet;
% cb=colorbar('TickLabelInterpreter', 'latex');
% % cb.Label.Interpreter='latex';
% set(gca,'XTick',0:0.1:1)
% axis([0.003 1 -inf inf])
% % 



% load fhn_03
% % subplot(1,2,2);
% d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
% er=pcolor(linspace(-0.09,0.09,501),linspace(0,2*pi,501),d);colormap jet;
% set(er,'EdgeColor','none')
% cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% % ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
% set(gca,'YTick',0:pi/4:2*pi);
% set(gca,'YTickLabel',ytl)
% axis square;ylabel('$\theta$');xlabel('$\Delta a$')


