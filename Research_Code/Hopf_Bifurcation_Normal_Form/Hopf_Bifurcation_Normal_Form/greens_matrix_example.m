close all
clear all
clc

tspan = 0:0.001:20;
x0=[1;0];
[t,x]=ode45(@(t,x) ms_system(t,x), tspan, x0);
n=2;
[~,locs]=findpeaks(x(:,1));
if length(locs)<1
    xg=x;
    tg=t;
    disp('Run Longer Simulation')
elseif length(locs)==1
    tg=t(1:locs(end))-t(1);
    xg=nan(length(tg),n);
    for ii=1:n
        xg(:,ii)=x(1:locs(end),ii);
    end
else
    tg=t(locs(end-1):locs(end))-t(locs(end-1));
    xg=nan(length(tg),n);
    for ii=1:n
        xg(:,ii)=x(locs(end-1):locs(end),ii);
    end
end

s.xg=xg;
s.tg=tg;

%%
p.tpert=0;
p.dt=0.001;
p.tau=1;

p.s0=[1 0 1 0];
[s.sensitivity_a, s.greens]=sensitivity_function(p,s);
[s.sensitivity]=variation_of_constants(p,s);

% plot(s.xg(:,1),s.xg(:,2),'r','linewidth',2);hold on;plot(s.xg(:,1),s.xg(:,2),'k','linewidth',2);plot(s.x_short(:,1),s.x_short(:,2),'g','linewidth',2)
% figure();subplot(2,1,1);plot(s.sensitivity(:,1),'linewidth',2);title('x');subplot(2,1,2);plot(s.sensitivity(:,2),'linewidth',2);title('y')





function xdot=ms_system(t,x)
xdot(1,1)=x(2);
xdot(2,1)=-x(1);
end 


