


a_vals=0:0.002:1;
tau_vals=round(linspace(0,s.tg(end)*0.2,length(a_vals)),3);
tau_vals(2)=0.002;
% tau_vals=round(s.tg(end)*0.2*(0:0.002:1),3);
% tau_vals=round((0:0.004:1)*s.tg(end),3);

final_vals=zeros(length(a_vals),length(tau_vals));
final_coords=cat(3,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals);
%[x y xp yp xr yr]
er=pcolor(a_vals,tau_vals,final_vals);colormap jet;
set(er,'EdgeColor','none')
pause(0.01)
kk=1;
total=length(1:length(tau_vals))*length(2:length(tau_vals));

for ii=1:length(a_vals)
    pausehere=1;
    for jj=length(tau_vals):-1:2
        kk=kk+1;
        if mod(kk,10)==0
            disp(round(100*kk/total,3))
        end
%         p.a=a_vals(ii);
%         p.dp=p.a(1)*0.5;
%         p.tau=tau_vals(jj);
        %%%
        %%%
        %%% note prc_a and prc are flipped here
%         p.x0=[sqrt(-p.a(1)/p.c(1)) 0];
%         numeric=0;
%         rr=sqrt(-p.a(1)/p.c(1));
%         T=round(2*pi/(p.b(1)+p.d(1)*(rr^2)),3);
%         t=0:2*pi/((T/p.dt)-1):2*pi;
%         s.xg=[(rr*cos(t))' (rr*sin(t))'];
%         s.tg=(0:p.dt:T)';
%         [s.prc, s.prc_a]=prc_function(p,s,numeric);
%         [s.irc,s.irc_a]=irc_function(p,s,numeric);
        %%%
        %%%
        %%%
        s.prc=s.prc_og;
        s.irc=s.irc_og;
        p.a=p.a(1);
        p.dp=a_vals(ii);
%         tau=round(s.tg(end)*0.2,3);
        p.tau=tau_vals(jj);
        p.tpert=0;
%         p.tpert=tau_vals(jj);
        p.tpert_idx=round(p.tpert/p.dt,5)+1;
        p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;
        if p.tau_idx>(s.tg(end)/p.dt)
            s.prc=[s.prc;[s.prc(:,1)+s.prc(end,1) s.prc(:,2:end)]];
            s.irc=[s.irc;[s.irc(:,1)+s.irc(end,1) s.irc(:,2:end)]];
        end
        p.perturbed=1;
        p.p_time=p.tpert:p.dt:p.tpert+p.tau;
        p.a=[p.a p.a+p.dp p.a];
        p.phi_span=2*pi*(p.tpert:p.dt:p.tpert+p.tau)/s.tg(end);
        [s.t_short,s.tp_short,s.x_short,s.xp_short]=hopf_ode(p);
        
        
        s.t_short=s.t_short(p.tpert_idx:end);
        s.tp_short=s.tp_short(p.tpert_idx:end);
        s.x_short=s.x_short(p.tpert_idx:end,:);
        s.xp_short=s.xp_short(p.tpert_idx:end,:);
        
        
        p.perturbed=0;
        %         [~,s.tp,~,s.xp]=hopf_ode(p);
        %% Sensitivity Function
        p.s0=[1 0 0 1];
        numeric_sens=0;
        [s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
        if numeric_sens==1
            [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
        else
            s.t_sensitivity=s.t_short;
            s.sensitivity=s.sensitivity_direct;
        end
        %% Coordinate Recovery
        s.xy_cr=cr_function(p,s,1);
        s.xy_recovered = s.x_short + s.xy_cr;
        %%
        final_vals(jj,ii)=sqrt((s.xy_recovered(end,1)-s.xp_short(end,1))^2+(s.xy_recovered(end,2)-s.xp_short(end,2))^2);
        final_coords(jj,ii,:)=[s.x_short(end,:) s.xp_short(end,:) s.xy_recovered(end,:)];
        set(er,'CData',final_vals);
        pause(0.0001)
    end
end
p.a=p.a(1);

% save hopf_04 final_coords

