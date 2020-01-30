% plot stuff

plot(s.xp_short(:,1),s.xp_short(:,2),'r','linewidth',2);title('trajectories')
hold on;plot(s.xg(:,1),s.xg(:,2),'k','linewidth',2);
plot(s.x_short(:,1),s.x_short(:,2),'g','linewidth',2)

figure();subplot(2,1,1);
plot(s.sensitivity(:,1),'linewidth',2);title('x');
title('sensitivity')
subplot(2,1,2);
plot(s.sensitivity(:,2),'linewidth',2);title('y')


figure();
for k=1:4
    subplot(2,2,k)
    title('Greens Matrix')
    hold on
    plot(s.greens(:,k),'k','linewidth',2)
    plot(s.greens_a(:,k),'b','linewidth',2)
    axis tight
end

figure()
jcob=nan(length(s.tp_short),4);
for k=1:length(s.tp_short)
x_p=s.xp_short(k,:)';
[D,~]=DF_nhopf(x_p,p);
jcob(k,:)=reshape(D,[1,4]);
end
for k=1:4;subplot(2,2,k);plot(jcob(:,k));title('Jacobian');end













% a=p.a(1);b=p.b;c=p.c;d=p.d;x=s.x_short(:,1);y=s.x_short(:,2);
% J=[a+c*(x.^2+y.^2)+2*x.*(c*x-d*y) -b-d*(x.^2+y.^2)+2*y.*(c*x-d*y) b+d*(x.^2+y.^2)+2*x.*(d*x+c*y) a+c*(x.^2+y.^2)+2*y.*(d*x+c*y)];
% for k=1:4;subplot(2,2,k);plot(J(:,k));title('Jacobian');end

%%%%%% in the variation of constants formula they integrate 
%%%%%% from 0 to t. Should this be p.tpert->p.tpert+p.tau?


% x02=[0.307999348196218 -0.131762671158053];
% [t,~,xsol,~]=hopf_ode(p);
% p.x0=x02;
% [t,~,xsol2,~]=hopf_ode(p);
% 
% hold on
% axis([min(xsol2(:,1)) max(xsol2(:,1))...
%     min(xsol2(:,2)) max(xsol2(:,2))])
% axis square; box on
% plot(s.xg(:,1),s.xg(:,2),'k','linewidth',2)
% line1=plot(xsol(1,1),xsol(1,2),'b.','markersize',20);
% line2=plot(xsol2(1,1),xsol2(1,2),'r.','markersize',20);
% pause(0.1)
% for k=1:length(t)
%     set(line1,'XData',xsol(k,1),'YData',xsol(k,2));
%     set(line2,'XData',xsol2(k,1),'YData',xsol2(k,2));
%     pause(0.001)
% end 
% 


% r=0.001:0.001:1.5*max(s.xg(:,1));
% theta=-a_phase(p,r,0);
% plot(s.xg(:,1),s.xg(:,2),'k','linewidth',2)
% hold on
% plot(r.*cos(theta),r.*sin(theta),'b','linewidth',2)

