function [irc_a, irc]=irc_function(p,s,numeric)
%% Analytic IRC
irc_a=nan;
%% Numeric IRC
irc=nan;
if numeric==1
    disp('Calculating IRC')
    V=s.V;
    scale_vector_idx=s.scale_vector_idx;
    xg=s.xg;
    tg=s.tg;
    
    xg_f=xg;
    
    ircx0=randi(100,1,3);
    tstep=tg(end):-p.dt:0;
    eps_2=10^-8;
    c0=ones(3,1)*ircx0;
    irc_pos_end=zeros(3);
    irc_neg_end=zeros(3);
    k=1;
    irc=[100 100 100;0 0 0;0 0 0];
    figure(100);
    subplot(3,1,1);line(1)=plot(0,0);subplot(3,1,2);line(2)=plot(0,0);subplot(3,1,3);line(3)=plot(0,0);
    pause(0.001)
    while sum(abs(irc(1,:)-irc(end,:))>10^-13)>0 && k<50
        k;
        if k>1
            c0=([1;1;1]*c0');
        end
        c0_pos=c0+eps_2*eye(3);
        c0_neg=c0-eps_2*eye(3);
        for n=1:3
            [~,irc_pos]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_pos(n,:));
            [~,irc_neg]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_neg(n,:));
            irc_pos_end(:,n)=irc_pos(end,:)';
            irc_neg_end(:,n)=irc_neg(end,:)';
        end
        c0_old=diag(c0);
        [~,irc]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_old);
%         irc=irc/(irc(1,:)*V(:,scale_vector_idx));
        J=eye(3)-((irc_pos_end-irc_neg_end))/(2*eps_2);
        c0=c0_old-J\(c0_old-irc(end,:)');
        k=k+1;
        set(line(1),'XData',1:length(irc(:,1)),'YData',irc(:,1))
        set(line(2),'XData',1:length(irc(:,1)),'YData',irc(:,2))
        set(line(3),'XData',1:length(irc(:,1)),'YData',irc(:,3))
        pause(0.001)
        abs(irc(1,:)-irc(end,:))
    end
    [t,irc]=(ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0));
    irc=irc/(irc(1,:)*V(:,scale_vector_idx));
    abs(irc(1,:)-irc(end,:));
    omega=2*pi/tg(end);
    t=omega*t;
    irc=flip([t irc]);
    set(line(1),'XData',1:length(irc(:,1)),'YData',irc(:,2))
    set(line(2),'XData',1:length(irc(:,1)),'YData',irc(:,3))
    set(line(3),'XData',1:length(irc(:,1)),'YData',irc(:,4))
    pause(0.001)
end




