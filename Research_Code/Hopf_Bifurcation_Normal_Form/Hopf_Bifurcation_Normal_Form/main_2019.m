close all
clear all
beep off
clc
set_default_plot
%% Parameters Setup
p=parameters(1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*25,3);
p.perturbed=0;  
pert=[1 0];
% TAU=20;
TAU=round(2*pi/((p.b(1)+p.d(1)*p.x0(1)^2))*0.2,3);
% TAU=1;
TPERT=0;
% p.tpert
% p.tau
%% 
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;

%% Do 2d stuff
% plots2d_approx
% figure();
% er=pcolor(a_vals,tau_vals,final_vals);colormap jet;
% set(er,'EdgeColor','none')

%% Perturbation 
p.tau=TAU;
p.tpert=TPERT;
p.tpert_idx=round(p.tpert/p.dt,5)+1;
p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;

if p.tau_idx>(s.tg(end)/p.dt)
    s.prc=[s.prc;[s.prc(:,1)+s.prc(end,1) s.prc(:,2:end)]];
    s.irc=[s.irc;[s.irc(:,1)+s.irc(end,1) s.irc(:,2:end)]];
end 

p.perturbed=1;
p.dp=p.a(1)*pert(1)+pert(2);
p.p_time=p.tpert:p.dt:p.tpert+p.tau;
p.a=[p.a p.a+p.dp p.a];
p.phi_span=2*pi*(p.tpert:p.dt:p.tpert+p.tau)/s.tg(end);
[s.t_short,s.tp_short,s.x_short,s.xp_short]=hopf_ode(p);


s.t_short=s.t_short(p.tpert_idx:end);
s.tp_short=s.tp_short(p.tpert_idx:end);
s.x_short=s.x_short(p.tpert_idx:end,:);
s.xp_short=s.xp_short(p.tpert_idx:end,:);


p.perturbed=0;
[~,s.tp,s.x,s.xp]=hopf_ode(p);
%% Sensitivity Function
p.s0=[1 0 0 1];
numeric_sens=0;
% disp('sens VOC 2')
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
% [s.sensitivity_direct, s.greens, s.greens_vector, s.greens_rs]=sensitivity_function_2(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
%     gt=flip(reshape(s.greens_vector',2,2,length(p.p_time)));
    for k=1:length(p.p_time)
        s.sensitivity(k,:)=(s.greens_rs(:,:,k)*s.sensitivity(k,:)')';
%         s.sensitivity(k,:)=(gt(:,:,k)*s.sensitivity(k,:)')';
    end
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end

% figure();
% plot(s.greens(:,1));hold on;plot(s.greens_vector(:,1))
% figure();
% plot(s.sensitivity);hold on;plot(s.sensitivity_direct)
%% Coordinate Recovery 
s.xy_cr=cr_function(p,s,0);
figure();hold on
plot(s.xg(:,1),s.xg(:,2))
plot(s.x_short(:,1),s.x_short(:,2))
plot(s.xp_short(:,1),s.xp_short(:,2))
s.xy_recovered = s.x_short + s.xy_cr;
plot(s.xy_recovered(:,1),s.xy_recovered(:,2))
% legend('x_g','x_s','x_p','x_r')
box on; 
mval=max(max(abs(s.xy_recovered)));
axis(1.3*[-mval mval -mval mval])
axis square
%%
% plot_stuff
% random_animation_script
%% Plot Phase Field 
[xm,ym]=meshgrid(-0.5:0.001:0.5,-0.5:0.001:0.5);
[m,n]=size(xm);
theta=nan(m,n);
for ii=1:m
    for jj=1:n
        theta(ii,jj)=a_phase(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2),...
            atan2(ym(ii,jj),xm(ii,jj)));
    end
end 
subplot(1,2,1)
pcolor(xm,ym,mod(theta,2*pi));shading flat
% colormap([jet;flipud(jet)])
cb=colorbar('TickLabelInterpreter', 'latex');set(cb,'Ticks',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
set(cb,'TickLabels',xtl)
colormap(jet);set(get(cb,'title'),'string','$\theta$','interpreter','latex');
axis square; axis tight;xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.125)
clear theta xm ym n m ii jj;hold on;plot(s.xg(:,1),s.xg(:,2),'w')
%% Plot Isostable Field 
% [xm,ym]=meshgrid([-0.5:0.001:-0.1 0.1:0.001:0.5],[-0.5:0.001:-0.1 0.1:0.001:0.5]);
[xm,ym]=meshgrid(-0.5:0.001:0.5,-0.5:0.001:0.5);
[m,n]=size(xm);
theta=nan(m,n);
for ii=1:m
    for jj=1:n
        if sqrt(xm(ii,jj)^2+ym(ii,jj)^2)>0.175
            theta(ii,jj)=a_approach(p,sqrt(xm(ii,jj)^2+ym(ii,jj)^2));
        end
    end
end
subplot(1,2,2)
pcolor(xm,ym,((theta)));shading flat
colormap(jet)
cb=colorbar('TickLabelInterpreter', 'latex')
colormap(jet);set(get(cb,'title'),'string','$\psi$','interpreter','latex');
axis square; axis tight;xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.125)
clear theta xm ym n m ii jj;hold on;plot(s.xg(:,1),s.xg(:,2),'w')


% 
% load hopf_06
% % subplot(1,2,2);
% d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
% er=pcolor(0:0.002:1,(0:0.002:1)*(2*pi/5),d);colormap jet;
% set(er,'EdgeColor','none')
% cb=colorbar('TickLabelInterpreter', 'latex');
% % ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
% set(gca,'YTick',0:pi/10:2*pi/5);
% set(gca,'YTickLabel',ytl)
% axis square;ylabel('$\theta$');xlabel('$a$')
% c2=caxis;
% c3 = [min([c1 c2]), max([c1 c2])];
% caxis(c3);
% subplot(1,2,1);colorbar off



