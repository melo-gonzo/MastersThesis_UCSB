function [irc_a, irc]=irc_function(p,s,numeric)
%% Analytic IRC
phi=2*pi*(s.tg)./s.tg(end);
r=sqrt(-p.a/p.c);
irc_a=[(phi)...
    (-sqrt(1+p.d(1)^2/p.c(1)^2)*cos(phi))...
    (-sqrt(1+p.d(1)^2/p.c(1)^2)*sin(phi))];

%% Numeric IRC
irc=nan;
if numeric==1
    disp('Calculating IRC')
    V=s.V;
    scale_vector_idx=s.scale_vector_idx;
    xg=s.xg;
    tg=s.tg;
    
    xg_f=xg;
    
    ircx0=rand(2,1)';
    tstep=tg(end):-p.dt:0;
    eps_2=10^-6;
    c0=ones(2,1)*ircx0;
    irc_pos_end=zeros(2);
    irc_neg_end=zeros(2);
    k=1;
    irc=[100 100;0 0];
    % sum(abs((irc(1,:)-irc(end,:)))>0.001)>0
    %
    while sum(abs((irc(1,:)-irc(end,:))>10^-5)) && k<50
        if k>1
            c0=([1;1]*c0');
        end
        c0_pos=c0+eps_2*eye(2);
        c0_neg=c0-eps_2*eye(2);
        for n=1:2
            [~,irc_pos]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_pos(n,:));
            [~,irc_neg]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_neg(n,:));
            irc_pos_end(:,n)=irc_pos(end,:)';
            irc_neg_end(:,n)=irc_neg(end,:)';
        end
        c0_old=diag(c0);
        [~,irc]=ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0_old);
        J=eye(2)-((irc_pos_end-irc_neg_end))/(2*eps_2);
        c0=c0_old-J\(c0_old-irc(end,:)');
        k=k+1;
    end
    [t,irc]=(ode45(@(t,irc_sol) irc_ode(t,irc_sol,xg_f,p,s),tstep,c0));
    irc=irc/(irc(1,:)*V(:,scale_vector_idx));
    
    abs(irc(1,:)-irc(end,:));
    omega=2*pi/tg(end);
    t=omega*t;
    irc=flip([t irc]);
end



