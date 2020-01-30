%% Figure 36

close all
clear all
beep off
clc
set_default_plot

ptrb=[0 5];
p=parameters(1,0.05,-70,3,50,5,-90,5,0,5);


p.x0=[-6.6507 0.24773 0.0017566];
p.dt=0.001;
p.tspan=0:p.dt:500;
p.perturbed=0;  

numeric=1;
[s.xg, s.tg, s.locs]=xg_function(p,numeric);
load prc 
load irc1
load irc2
s.prc=prc;s.irc=irc1;s.irc2=irc2;

s.prc_og=s.prc;
s.irc_og=s.irc;
s.irc_og2=s.irc2;
p.dp=p.ib(1)*ptrb(1)+ptrb(2);
plots2d_approx


figure();







h(1) = subplot(2,2,1);
er=pcolor((1:1680)*TAU/1680,linspace(0,1,length(tau_vals))*(2*pi),final_vals_theta);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};

set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \theta$','interpreter','latex');

h(2) = subplot(2,2,2);
er=pcolor((1:1680)*TAU/1680,linspace(0,1,length(tau_vals))*(2*pi),final_vals_psi);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \psi_1$','interpreter','latex');


h(3) = subplot(2,2,3);
er=pcolor((1:1680)*TAU/1680,linspace(0,1,length(tau_vals))*(2*pi),final_vals_psi2);colormap jet;
set(er,'EdgeColor','none')
cb=colorbar('TickLabelInterpreter', 'latex');
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
ytl={'$0$','$\frac{\pi}{2}$','$\pi$','$\frac{3\pi}{2}$','$2\pi$'};
% ytl={'$0$','$\frac{\pi}{10}$','$\frac{\pi}{5}$','$\frac{3\pi}{10}$','$\frac{2\pi}{5}$'};
set(gca,'YTick',0:pi/2:2*pi);
set(gca,'YTickLabel',ytl)
axis square;ylabel('$\theta_0$');xlabel('$\tau$');rotate_y_label(0,-0.2)
set(get(cb,'title'),'string','$\Delta \psi_2$','interpreter','latex');

pos = get(h,'Position');
new = mean(cellfun(@(v)v(1),pos(1:2)));
set(h(3),'Position',[new,pos{end}(2:end)])





%%

% h(1) = subplot(2,2,1);
% h(2) = subplot(2,2,2);
% h(3) = subplot(2,2,3);
% pos = get(h,'Position');
% new = mean(cellfun(@(v)v(1),pos(1:2)));
% set(h(3),'Position',[new,pos{end}(2:end)])


