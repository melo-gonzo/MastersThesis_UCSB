%% Figure 1, 
close all
clear all
beep off
clc
set_default_plot
p=parameters(.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=0:p.dt:100;
p.perturbed=0;  
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
% p.tau=round(s.tg(end)*0.2,3);
p.tau=round(s.tg(end)*0.2,3);
p.tpert=round(s.tg(end)*1.6,3);
p.tpert_idx=round(p.tpert/p.dt,5)+1;
p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;
p.dp=p.a(1)*0+1;
p.p_time=p.tpert:p.dt:p.tpert+p.tau;
p.a=[p.a p.a+p.dp p.a];
[s.t_short,s.tp_short,s.x_short,s.xp_short]=hopf_ode(p);
p.perturbed=0;
[~,s.tp,s.x,s.xp]=hopf_ode(p);
hold on
plot(s.tp,s.x(:,1))
plot(s.tp,s.xp(:,1),'b-.')
ps=s.tp(p.tpert_idx);
pf=s.tp(p.tau_idx);
p1=line([ps ps],[min(s.xp(:,1)) max(s.xp(:,1))]);
p2=line([pf pf],[min(s.xp(:,1)) max(s.xp(:,1))]);
set([p1 p2],'visible','off')
gray=0.7;
patch([ps pf pf ps],[min(s.xp(:,1)) min(s.xp(:,1))...
    max(s.xp(:,1)) max(s.xp(:,1))],[gray gray gray])
set(gca,'children',flipud(get(gca,'children')))
axis([0 s.tg(end)*5 -inf inf])
xlabel('$t$');ylabel('$x$');rotate_y_label(0,-0.075)
get(gcf,'Position')
set(gcf,'Position',[ans(1:3) ans(4)/2])

% figure();
% p.tpert_idx=round(p.tpert/p.dt,5)+1;
% p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;
% s.t_short=s.t_short(p.tpert_idx:end);
% s.tp_short=s.tp_short(p.tpert_idx:end);
% s.x_short=s.x_short(p.tpert_idx:end,:);
% s.xp_short=s.xp_short(p.tpert_idx:end,:);
% plot(s.xg(:,1),s.xg(:,2))
% plot(s.xp_short(:,1),s.xp_short(:,2),'b-.')
% xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.075)

%% Figure 9
close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=0:p.dt:100;
p.perturbed=0;  
pert=[0 0.1];
TAU=1;
TPERT=0;
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
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
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
    for k=1:length(p.p_time)
        s.sensitivity(k,:)=(s.greens_rs(:,:,k)*s.sensitivity(k,:)')';
    end
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end
%%% Go in here and pause to plot.
s.xy_cr=cr_function(p,s,0);

%% Figure 6 and 7
close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=0:p.dt:100;
p.perturbed=0;  
pert=[0 0.1];
TAU=1;
TPERT=0;
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);



figure();
subplot(1,2,1);plot(s.prc(:,1),s.prc(:,2))
ylabel('$\frac{\partial \theta}{\partial x}$','FontSize',12);
xlabel('$\theta$')
rotate_y_label(0,-0.15)
set(gca,'XTick',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
set(gca,'XTickLabel',xtl)
axis([0 2*pi -inf inf])
subplot(1,2,2);plot(s.prc(:,1),s.prc(:,3))
ylabel('$\frac{\partial \theta}{\partial y}$','FontSize',12);
rotate_y_label(0,-0.15)
xlabel('$\theta$')
set(gca,'XTick',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
set(gca,'XTickLabel',xtl)
axis([0 2*pi -inf inf])
half_figure(0.5)

figure();
s.irc(:,2:end)=-s.irc(:,2:end);
subplot(1,2,1);plot(s.irc(:,1),s.irc(:,2))
ylabel('$\frac{\partial \psi}{\partial x}$','FontSize',12);
xlabel('$\theta$')
rotate_y_label(0,-0.15)
set(gca,'XTick',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
set(gca,'XTickLabel',xtl)
axis([0 2*pi -inf inf])
subplot(1,2,2);plot(s.irc(:,1),s.irc(:,3))
ylabel('$\frac{\partial \psi}{\partial y}$','FontSize',12);
rotate_y_label(0,-0.15)
xlabel('$\theta$')
set(gca,'XTick',0:pi/4:2*pi)
xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
set(gca,'XTickLabel',xtl)
axis([0 2*pi -inf inf])
half_figure(0.5)

%% Figure 15
close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=0:p.dt:100;
p.perturbed=0;  
pert=[0 0];
TAU=round(2*pi/((p.b(1)+p.d(1)*p.x0(1)^2))*3,3);
TPERT=0;
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
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
s.perturbed=0;
[~,s.tp,s.x,s.xp]=hopf_ode(p);
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
    for k=1:length(p.p_time)
        s.sensitivity(k,:)=(s.greens_rs(:,:,k)*s.sensitivity(k,:)')';
    end
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end
figure();
subplot(1,2,1);plot(s.t_sensitivity,s.sensitivity(:,1))
ylabel('$\frac{\partial x}{\partial a}$','FontSize',12);
xlabel('$\tau$')
rotate_y_label(0,-0.15)
axis tight; 
% axis square
subplot(1,2,2);plot(s.t_sensitivity,s.sensitivity(:,2))
ylabel('$\frac{\partial y}{\partial a}$','FontSize',12);
rotate_y_label(0,-0.15)
xlabel('$\tau$')
half_figure(0.5)
axis tight; 
% axis square 
% set(gcf,'Position',[940 505 440 420])
% half_figure(0.5);

%% Figure 10
close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*25,3);
p.perturbed=0;  
pert=[0 0.3];
TAU=round(2*pi/((p.b(1)+p.d(1)*p.x0(1)^2))*0.2,3);
TAU=0.8;
TPERT=0;
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
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

figure();hold on
plot(s.xg(:,1),s.xg(:,2))
plot(s.x_short(:,1),s.x_short(:,2),'g')
plot(s.x_short(end,1),s.x_short(end,2),'g.','markersize',20)
plot(s.xp_short(:,1),s.xp_short(:,2),'b-.')
plot(s.xp_short(end,1),s.xp_short(end,2),'b.','markersize',20)
% axis square
axis([0.15 0.35 -0.01 0.35])
xlabel('$x$')
ylabel('$y$')
win=get(gca,'Position');
x1=[win(1)+((s.x_short(end,1)-0.15)/(.35-.15))*(win(3)) win(1)+((s.xp_short(end,1)-0.15)/(.35-.15))*(win(3))];
x2=[win(2)+((s.x_short(end,2)+0.01)/(.35+0.01))*(win(4)) win(2)+((s.x_short(end,2)+0.01)/(.35+0.01))*(win(4))];
a = annotation('textarrow',x1,x2)
x1=[win(1)+((s.x_short(end,1)-0.15)/(.35-.15))*(win(3)) win(1)+((s.x_short(end,1)-0.15)/(.35-.15))*(win(3))];
x2=[win(2)+((s.x_short(end,2)+0.01)/(.35+0.01))*(win(4)) win(2)+((s.xp_short(end,2)+0.01)/(.35+0.01))*(win(4))];
a = annotation('textarrow',x1,x2)
text(0.22,0.2325,'$\Delta x$')
text(0.1875,0.2725,'$\Delta y$')
set(gcf,'Position',[940 505 440 420])

%% Figure 11
close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*25,3);
p.perturbed=0;  
pert=[0 0.3];
TAU=0.8;
TPERT=0;
numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
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
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
if numeric_sens==1
    [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
    for k=1:length(p.p_time)
        s.sensitivity(k,:)=(s.greens_rs(:,:,k)*s.sensitivity(k,:)')';
    end
else
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
end

s.xy_cr=cr_function(p,s,0);
figure();hold on
plot(s.xg(:,1),s.xg(:,2))
% plot(s.x_short(:,1),s.x_short(:,2))
plot(s.xp_short(:,1),s.xp_short(:,2),'b-.')
s.xy_recovered = s.x_short + s.xy_cr;
plot(s.xy_recovered(:,1),s.xy_recovered(:,2),'r:')
axis([0.19 0.55 -0.01 0.4])
xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.1);
box on; 
axis square

axes('Position',[0.48 0.5 0.395 0.4])
hold on;plot(s.xp_short(:,1),s.xp_short(:,2),'b-.');plot(s.xy_recovered(:,1),s.xy_recovered(:,2),'r:')
axis([0.241 0.246 0.299 0.31])
xlabel('x');ylabel('y');
text(0.24275,0.3035,'$d=\sqrt{\Delta \tilde{x}^2 + \Delta \tilde{y}^2}$');
% axis square

win=get(gca,'Position')
xlim=get(gca,'XLim');
ylim=get(gca,'YLim');

x1=[win(1)+((s.xy_recovered(end,1)-xlim(1))/(xlim(2)-xlim(1)))*(win(3))...
    win(1)+((s.xp_short(end,1)-xlim(1))/(xlim(2)-xlim(1)))*(win(3))];

x2=[win(2)+((s.xy_recovered(end,2)-ylim(1))/(ylim(2)-ylim(1)))*(win(4))...
    win(2)+((s.xp_short(end,2)-ylim(1))/(ylim(2)-ylim(1)))*(win(4))];
a = annotation('textarrow',x1,x2)
xlabel('$x$');ylabel('$y$');rotate_y_label(0,-0.225);
set(gcf,'Position',[940 505 440 420])

%% Figure 12
close all
clear all
clc
load hopf_01
d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
er=pcolor(0:0.004:1,(0:0.004:1)*(2*pi),d);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
set(gca,'YTick',0:pi/4:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta$');xlabel('$\Delta a$');rotate_y_label(0,-0.1)

%% Figure 13
% close all;clear all;clc
load hopf_03
subplot(1,2,2);
d=sqrt((final_coords(:,:,3)-final_coords(:,:,5)).^2 +(final_coords(:,:,4)-final_coords(:,:,6)).^2);
er=pcolor(0:0.002:1,(0:0.002:1)*(2*pi/5),d);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/10:2*pi/5);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta$');xlabel('$\Delta a$');rotate_y_label(0,-0.175)
% c2=caxis;
% c3 = [min([c1 c2]), max([c1 c2])];
% caxis(c3);
% subplot(1,2,1);colorbar off

%% Figure 34

close all
clear all
beep off
clc
set_default_plot
p=parameters(0.1,1,-1,1);
p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
p.dt=0.001;
p.tspan=round(0:p.dt:5.712*20,3);
p.perturbed=0;  
pert=[1 0];

TAU=round(2*pi/((p.b(1)+p.d(1)*p.x0(1)^2))*0.2,3);
% TAU=1;
TPERT=0;

numeric=1;
[s.xg, s.tg]=xg_function(p,numeric);
[s.prc_a, s.prc]=prc_function(p,s,numeric);
[s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
[s.irc_a,s.irc]=irc_function(p,s,numeric);
s.prc_og=s.prc;
s.irc_og=s.irc;
p.dp=p.a(1)*pert(1)+pert(2);
plots2d_approx

subplot(1,2,1);
er=pcolor((1:1143)*TAU/1143,linspace(0,1,length(tau_vals))*(2*pi),final_vals_theta);colormap jet;
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
er=pcolor((1:1143)*TAU/1143,linspace(0,1,length(tau_vals))*(2*pi),final_vals_psi);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \psi$','interpreter','latex');



