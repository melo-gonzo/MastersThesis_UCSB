function [prc_a, prc] = prc_function(p,s,numeric)
xg=s.xg;
tg=s.tg;
%% Analytic PRC
prc_a=nan;
%% Numeric PRC
prc=nan;
if numeric==1
    disp('Calculating PRC')
    options=odeset('RelTol',3e-14,'AbsTol',1e-17*ones(length(p.x0),1));
    % Solve for periodic solution backwards in time using newton method
    % xg_f=flip(xg);
    xg_f=xg;
    prcx0=randi(10,1,3);
    % prcx0=[4 2];
    tstep=tg(end):-p.dt:0;
    eps_2=10^-6;
    c0=ones(3,1)*prcx0;
    prc_pos_end=zeros(3);
    prc_neg_end=zeros(3);
    k=1;
    prc=[100 100 100;0 0 0;0 0 0];
    % sum(abs((prc(1,:)-prc(end,:))./prc(1,:))>0.001)>0
    figure(100);
    subplot(3,1,1);line(1)=plot(0,0);subplot(3,1,2);line(2)=plot(0,0);subplot(3,1,3);line(3)=plot(0,0);
    pause(0.001)
    while sum(abs(prc(1,:)-prc(end,:))>10^-15)>0 && k<50
        k;
        if k>1
            c0=([1;1;1]*c0');
        end
        c0_pos=c0+eps_2*eye(3);
        c0_neg=c0-eps_2*eye(3);
        for n=1:3
            [~,prc_pos]=ode45(@(t,prc_sol) prc_ode(t,prc_sol,xg_f,p,s),tstep,c0_pos(n,:));
            [~,prc_neg]=ode45(@(t,prc_sol) prc_ode(t,prc_sol,xg_f,p,s),tstep,c0_neg(n,:));
            prc_pos_end(:,n)=prc_pos(end,:)';
            prc_neg_end(:,n)=prc_neg(end,:)';
        end
        c0_old=diag(c0);
        [~,prc]=ode45(@(t,prc_sol) prc_ode(t,prc_sol,xg_f,p,s),tstep,c0_old);
        J=eye(3)-((prc_pos_end-prc_neg_end))/(2*eps_2);
        c0=c0_old-J\(c0_old-prc(end,:)');
        k=k+1;
        set(line(1),'XData',1:length(prc(:,1)),'YData',prc(:,1))
        set(line(2),'XData',1:length(prc(:,1)),'YData',prc(:,2))
        set(line(3),'XData',1:length(prc(:,1)),'YData',prc(:,3))
        pause(0.001);
    end
    [t,prc]=(ode45(@(t,prc_sol) prc_ode(t,prc_sol,xg_f,p,s),tstep,c0));
    t=flip(t);prc=flip(prc);
    omega=2*pi/tg(end);
    prc=prc/(xg(1,:)*prc(1,:)');
    t=omega*t;
    prc=[t prc];
    set(line(1),'XData',1:length(prc(:,1)),'YData',prc(:,2))
    set(line(2),'XData',1:length(prc(:,1)),'YData',prc(:,3))
    set(line(3),'XData',1:length(prc(:,1)),'YData',prc(:,4))
    pause(0.001);
    
end




