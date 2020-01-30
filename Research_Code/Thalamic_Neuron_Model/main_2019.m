close all
clear all
beep off
clc
set_default_plot
%% Parameters Setup
ptrb=[0 15];
% cm,gl,el,gna,ena,gk,ek,gt,et,ib
p=parameters(1,0.05,-70,3,50,5,-90,5,0,5);
TPERT=2;

p.x0=[-6.6507 0.24773 0.0017566];
p.dt=0.001;
p.tspan=0:p.dt:500;
p.perturbed=0;  
% [s.t,~,s.x,~]=thalamic_ode(p);

%% 
numeric=1;
[s.xg, s.tg, s.locs]=xg_function(p,numeric);
% xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
% subplot(1,3,1);plot(s.tg,s.xg(:,1));set(gca,'XTick',0:s.tg(end)/8:s.tg(end));set(gca,'XTickLabel',xtl);xlabel('$\theta$');ylabel('$v$');axis tight;axis square;rotate_y_label(0,-0.225)
% subplot(1,3,2);plot(s.tg,s.xg(:,2));set(gca,'XTick',0:s.tg(end)/8:s.tg(end));set(gca,'XTickLabel',xtl);xlabel('$\theta$');ylabel('$h$');axis tight;axis square;rotate_y_label(0,-0.225)
% subplot(1,3,3);plot(s.tg,s.xg(:,3));set(gca,'XTick',0:s.tg(end)/8:s.tg(end));set(gca,'XTickLabel',xtl);xlabel('$\theta$');ylabel('$r$');axis tight;axis square;rotate_y_label(0,-0.225)
% half_figure(0.5)
% [s.prc_a, s.prc]=prc_function(p,s,numeric);
% [s.kappa, s.V, s.scale_vector_idx]=floquet_function(p,s,numeric);
% s.scale_vector_idx=2;
% [s.irc_a,s.irc]=irc_function(p,s,numeric);
% s.scale_vector_idx=1;
% [s.irc_a2,s.irc2]=irc_function(p,s,numeric);

load prc 
load irc1
load irc2
s.prc=prc;s.irc=irc1;s.irc2=irc2;


% close all
% movetitle=-0.2
% xtl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'}
% subplot(3,2,1);plot(s.irc(:,1),s.irc(:,2));ylabel('$\frac{\partial \psi_1}{\partial v}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');
% axis tight;rotate_y_label(0,movetitle)
% subplot(3,2,3);plot(s.irc(:,1),s.irc(:,3));ylabel('$\frac{\partial \psi_1}{\partial h}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');
% axis tight;rotate_y_label(0,movetitle)
% subplot(3,2,5);plot(s.irc(:,1),s.irc(:,4));ylabel('$\frac{\partial \psi_1}{\partial r}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');rotate_y_label(0,movetitle)
% axis([-inf 2*pi -inf inf]);axis tight
% 
% subplot(3,2,2);plot(s.irc2(:,1),s.irc2(:,2));ylabel('$\frac{\partial \psi_2}{\partial v}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');
% axis tight;rotate_y_label(0,movetitle)
% subplot(3,2,4);plot(s.irc2(:,1),s.irc2(:,3));ylabel('$\frac{\partial \psi_2}{\partial h}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');
% axis tight;rotate_y_label(0,movetitle)
% subplot(3,2,6);plot(s.irc2(:,1),s.irc2(:,4));ylabel('$\frac{\partial \psi_2}{\partial r}$','FontSize',12);
% set(gca,'XTick',0:pi/4:2*pi);set(gca,'XTickLabel',xtl);xlabel('$\theta$');rotate_y_label(0,movetitle)
% axis([-inf 2*pi -inf inf]);axis tight



s.prc_og=s.prc;
s.irc_og=s.irc;
s.irc_og2=s.irc2;

% figure();
% subplot(3,1,1);title('PRC');plot([s.prc(:,2);s.prc(:,2)]);ylabel('v');subplot(3,1,2);plot([s.prc(:,3);s.prc(:,3)]);ylabel('h');subplot(3,1,3);plot([s.prc(:,4);s.prc(:,4)]);ylabel('r')
% figure()
% subplot(3,1,1);title('IRC_1');plot([s.irc(:,2);s.irc(:,2)]);ylabel('v');subplot(3,1,2);plot([s.irc(:,3);s.irc(:,3)]);ylabel('h');subplot(3,1,3);plot([s.irc(:,4);s.irc(:,4)]);ylabel('r')
% figure()
% subplot(3,1,1);title('IRC_2');plot([s.irc2(:,2);s.irc2(:,2)]);ylabel('v');subplot(3,1,2);plot([s.irc2(:,3);s.irc2(:,3)]);ylabel('h');subplot(3,1,3);plot([s.irc2(:,4);s.irc2(:,4)]);ylabel('r')
%% Do 2d stuff
plots2d
figure();
er=pcolor(a_vals,tau_vals,abs(final_vals));colormap jet;
set(er,'EdgeColor','none')

%% Perturbation 
p.x0=s.xg(1,:);
p.tau=round(s.tg(end)*0.2,3);

% p.tpert=round(s.tg(end)*TPERT/(2*pi),3);
% p.tpert=round(s.tg(end)*4.86/(2*pi),3);
p.tpert=TPERT;

p.tpert_idx=round(p.tpert/p.dt,5)+1;
p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;

if p.tau_idx>(s.tg(end)/p.dt)
    s.prc=[s.prc;[s.prc(:,1)+s.prc(end,1) s.prc(:,2:end)]];
    s.irc=[s.irc;[s.irc(:,1)+s.irc(end,1) s.irc(:,2:end)]];
    s.irc2=[s.irc2;[s.irc2(:,1)+s.irc2(end,1) s.irc2(:,2:end)]];
end 

p.perturbed=1;
p.dp=p.ib(1)*ptrb(1)+ptrb(2);
p.p_time=p.tpert:p.dt:p.tpert+p.tau;
p.ib=[p.ib p.ib+p.dp p.ib];
p.phi_span=2*pi*(p.tpert:p.dt:p.tpert+p.tau)/s.tg(end);
[s.t_short,s.tp_short,s.x_short,s.xp_short]=thalamic_ode(p);
s.t_short=s.t_short(p.tpert_idx:end);
s.tp_short=s.tp_short(p.tpert_idx:end);
s.x_short=s.x_short(p.tpert_idx:end,:);
s.xp_short=s.xp_short(p.tpert_idx:end,:);
p.perturbed=0;
[~,s.tp,s.x,s.xp]=thalamic_ode(p);
s.xg_p=s.xp(s.locs(1):s.locs(2),:);

%% Sensitivity Function
p.s0=[1 0 0 1];
numeric_sens=0;
[s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
s.t_sensitivity=s.t_short;
s.sensitivity=s.sensitivity_direct;


% figure();movetitle=-0.075;
% subplot(3,1,1);hold on;plot(s.t_sensitivity,s.sensitivity(:,1));ylabel('$\frac{\partial v}{\partial I_b}$','FontSize',12);axis tight;rotate_y_label(0,movetitle)
% subplot(3,1,2);hold on;plot(s.t_sensitivity,s.sensitivity(:,2));ylabel('$\frac{\partial h}{\partial I_b}$','FontSize',12);axis tight;rotate_y_label(0,movetitle)
% subplot(3,1,3);hold on;plot(s.t_sensitivity,s.sensitivity(:,3));ylabel('$\frac{\partial r}{\partial I_b}$','FontSize',12);axis tight;xlabel('$\tau$');rotate_y_label(0,movetitle)


%% Coordinate Recovery 
s.xy_cr=cr_function(p,s,0);
s.xy_recovered = s.x_short + s.xy_cr;
%% Plot Coordinate Recovery Plots
figure();hold on
plot3(s.xg(:,1),s.xg(:,2),s.xg(:,3))
plot3(s.x_short(:,1),s.x_short(:,2),s.x_short(:,3))
plot3(s.xp_short(:,1),s.xp_short(:,2),s.xp_short(:,3),'b-.')
%% Coordinate Recovery Plot
plot3(s.xy_recovered(:,1),s.xy_recovered(:,2),s.xy_recovered(:,3))
view([-15 15])

figure();
for k=1:3
subplot(3,1,k);
hold on;
plot(s.tg,s.xg(:,k))
plot(s.t_short,s.x_short(:,k))
plot(s.t_short,s.xp_short(:,k))
plot(s.t_short,s.xy_recovered(:,k))
plot(s.tg,s.xg_p(:,k),'b-.')
axis tight
end

% xlim=[ -66.11      -6.6507];% figure()
% load thal_01
% d=final_coords(:,:,end);
% a_vals=0:0.003:1;
% d=pcolor(0:0.05:25,linspace(0,2*pi,length(0:0.05:25)),d)
% set(d,'EdgeColor','none');colormap jet;
% cb=colorbar('TickLabelInterpreter', 'latex');
% xlabel('$\Delta I_b$')
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% set(gca,'YTick',0:pi/4:2*pi)
% set(gca,'YTickLabel',ytl)
% ylabel('$\theta$')
% rotate_y_label(0,-0.075)

% ylim=[0.06549      0.80621];
% zlim=[0.0017012    0.0019149];
% set(gca,'Ylim',ylim)
% set(gca,'XLim',xlim)
%% YE
% figure()
% load thal_01
% % d=final_coords(:,:,end);
% a_vals=0:0.003:1;
% d=sqrt((final_coords(:,:,7)-final_coords(:,:,4)).^2+(final_coords(:,:,8)-final_coords(:,:,5)).^2+(final_coords(:,:,6)-final_coords(:,:,9)).^2);
% d=pcolor(0:0.05:25,linspace(0,2*pi,length(0:0.05:25)),d)
% set(d,'EdgeColor','none');colormap jet;
% cb=colorbar('TickLabelInterpreter', 'latex');
% xlabel('$\Delta I_b$')
% ytl={'$0$','$\frac{\pi}{4}$','$\frac{\pi}{2}$','$\frac{3\pi}{4}$','$\pi$','$\frac{5\pi}{4}$','$\frac{3\pi}{2}$','$\frac{7\pi}{4}$','$2\pi$'};
% set(gca,'YTick',0:pi/4:2*pi)
% set(gca,'YTickLabel',ytl)
% ylabel('$\theta$')
% rotate_y_label(0,-0.075)
% axis square
% % 






