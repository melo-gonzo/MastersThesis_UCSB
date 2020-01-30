%% Figure 16
close all
clear all
beep off
clc
set_default_plot
epss=0.2;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
subplot(1,3,1);
plot(s.tg(:,1),s.xg(:,1));
set(gca,'XTick',0:s.tg(end)/8:s.tg(end));
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
set(gca,'XTickLabel',xtl)
ylabel('$v$');rotate_y_label(0,-0.2)
xlabel('$\theta$');axis tight; axis square
subplot(1,3,2);
plot(s.tg(:,1),s.xg(:,2));
set(gca,'XTick',0:s.tg(end)/8:s.tg(end));
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
set(gca,'XTickLabel',xtl)
ylabel('$w$');rotate_y_label(0,-0.2)
xlabel('$\theta$');axis tight;axis square
subplot(1,3,3);
plot(s.xg(:,1),s.xg(:,2))
ylabel('$v$');rotate_y_label(0,-0.2)
xlabel('$w$');axis tight;axis square
half_figure(0.5)
aaaa=1;
%% Figure 17
close all
clear all
beep off
clc
set_default_plot

epss=0.80;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}

movey=-0.2;
subplot(3,2,1)
plot(s.prc(:,1),s.prc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square

subplot(3,2,2)
plot(s.prc(:,1),s.prc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square





epss=0.50;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);

subplot(3,2,3)
plot(s.prc(:,1),s.prc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])

subplot(3,2,4)
plot(s.prc(:,1),s.prc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])
% set(gcf,'Position',[940 505 440 420])










epss=0.20;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);

subplot(3,2,5)
plot(s.prc(:,1),s.prc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])

subplot(3,2,6)
plot(s.prc(:,1),s.prc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \theta}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])


pspsp=[956   173   463   632]

set(gcf,'Position',pspsp)
% 948        270.5          456        564.5
%% Figure 19
close all
clear all
beep off
clc
set_default_plot
epss=0.80;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}

movey=-0.2;
subplot(3,2,1)
plot(s.irc(:,1),s.irc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square

subplot(3,2,2)
plot(s.irc(:,1),s.irc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square






epss=0.50;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
subplot(3,2,3)
plot(s.irc(:,1),s.irc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])

subplot(3,2,4)
plot(s.irc(:,1),s.irc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])




















epss=0.20;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
subplot(3,2,5)
plot(s.irc(:,1),s.irc(:,2));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial v}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])

subplot(3,2,6)
plot(s.irc(:,1),s.irc(:,3));
set(gca,'XTick',0:pi/4:2*pi);
set(gca,'XTickLabel',xtl)
ylabel('$\frac{\partial \psi}{\partial w}$','FontSize',12);rotate_y_label(0,movey);
xlabel('$\theta$'); axis tight; axis square;axis([0 2*pi -inf inf])



pspsp=[956   173   463   632];
set(gcf,'Position',pspsp)


%% Figure 21
% close all
clear all
beep off
clc
set_default_plot
epss=0.20;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
p.x0=s.xg(1,:);
p.tau=round(s.tg(end)*3,3);
p.tpert=0;
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
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end

subplot(1,2,1);hold on
plot(s.t_sensitivity,s.sensitivity(:,1));
ylabel('$\frac{\partial v}{\partial a}$','Fontsize',12,'interpreter','latex');rotate_y_label(0,-0.175)
axis tight;xlabel('$\tau$');
% axis square
subplot(1,2,2);hold on
plot(s.t_sensitivity,s.sensitivity(:,2));
ylabel('$\frac{\partial w}{\partial a}$','Fontsize',12);rotate_y_label(0,-0.175)
axis tight;xlabel('$\tau$');
% axis square;
% set(gcf,'Position',[940 505 440 420])
half_figure(0.5);

%% Figure 22

close all
clear all
beep off
clc
set_default_plot
epss=0.50;
ptrb=[0 -0.15];
p=parameters(0,0.8,0,0,epss);
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

subplot(1,2,1);
plot(xg1(:,1),xg1(:,2))
hold on;plot(xg2(:,1),xg2(:,2),'b-.')
ylabel('$v$');rotate_y_label(0,-0.175); 
% axis square;
xlabel('$w$')

ptrb=[0 0.15];
p=parameters(0,0.8,0,0,epss);
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

subplot(1,2,2);
plot(xg1(:,1),xg1(:,2))
hold on;plot(xg2(:,1),xg2(:,2),'b-.')
ylabel('$v$');rotate_y_label(0,-0.175); 
% axis square;
xlabel('$w$')
half_figure(0.5)
% set(gcf,'Position',[940 505 440 420])

%% Figure 24


load fhn_01
% subplot(1,2,2);
d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
er=pcolor(linspace(-0.23,0.23,501),linspace(0,2*pi,501),d);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/4:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta$');xlabel('$\Delta a$'); rotate_y_label(0,-0.1)

%% Figure 32

close all
clear all
beep off
clc
set_default_plot
epss=0.50;
ptrb=[0 0.09];
p=parameters(0,0.8,0,0,epss);
p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
p.x0=s.xg(1,:);
p.tau=round(s.tg(end)*0.2,3);
% p.tpert=0;
p.tpert=round(s.tg(end)*pi/2/(2*pi),3);
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
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end
s.xy_cr=cr_function(p,s,0);
s.xy_recovered = s.x_short + s.xy_cr;
figure(3);hold on
plot(s.xg(:,1),s.xg(:,2))
plot(s.x_short(:,1),s.x_short(:,2))
plot(s.xp_short(:,1),s.xp_short(:,2))
plot(s.xy_recovered(:,1),s.xy_recovered(:,2))

%% Figure 35

close all
clear all
beep off
clc
set_default_plot
epss=0.50;
ptrb=[0 0.1];
p=parameters(0,0.8,0,0,epss);
TPERT=0;

p.x0=[-1 0];
p.dt=0.001;
p.tspan=0:p.dt:200;
p.perturbed=0;  
[s.t,~,s.x,~]=fhn_ode(p);

numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
p.dp=p.a(1)*ptrb(1)+ptrb(2);

plots2d_approx
figure();

subplot(1,2,1);
er=pcolor((1:2371)*TAU/2371,linspace(0,1,length(tau_vals))*(2*pi),final_vals_theta);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};

set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \theta$','interpreter','latex');

subplot(1,2,2);
er=pcolor((1:2371)*TAU/2371,linspace(0,1,length(tau_vals))*(2*pi),final_vals_psi);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \psi$','interpreter','latex');

