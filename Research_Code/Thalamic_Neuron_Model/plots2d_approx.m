


% a=round(linspace(0,s.tg(end)*0.2,length(a_vals)),3);
TAU=round(s.tg(end)*0.2,3);
% a(2)=0.002;R\
tau_vals=round((0:0.004:1)*s.tg(end),3);

final_vals_theta=zeros(length(tau_vals),1680);
final_vals_psi=zeros(length(tau_vals),1680);
final_vals_psi_2=zeros(length(tau_vals),1680);
% final_coords=cat(3,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals);
%[x y xp yp xr yr]
er=pcolor((1:1680),tau_vals,final_vals_theta);colormap jet;
set(er,'EdgeColor','none')
pause(0.01)
kk=1;
total=length(1:length(tau_vals))*length(2:length(tau_vals));

for jj=length(tau_vals):-1:1
    kk=kk+1;
    if mod(kk,10)==0
        disp(round(100*kk/total,3))
    end
    
    s.prc=s.prc_og;
    s.irc=s.irc_og;
    s.irc2=s.irc_og2;
    p.ib=p.ib(1);
    p.dp=p.dp(1);
    p.tau=TAU;
    p.tpert=tau_vals(jj);
    %         p.tpert=tau_vals(jj);
    p.tpert_idx=round(p.tpert/p.dt,5)+1;
    p.tau_idx=round((p.tpert+p.tau)/p.dt,5)+1;
    if p.tau_idx>(s.tg(end)/p.dt)
        s.prc=[s.prc;[s.prc(:,1)+s.prc(end,1) s.prc(:,2:end)]];
        s.irc=[s.irc;[s.irc(:,1)+s.irc(end,1) s.irc(:,2:end)]];
        s.irc2=[s.irc2;[s.irc2(:,1)+s.irc2(end,1) s.irc2(:,2:end)]];
    end
    p.perturbed=1;
    p.p_time=p.tpert:p.dt:p.tpert+p.tau;
    p.ib=[p.ib p.ib+p.dp p.ib];
    p.phi_span=2*pi*(p.tpert:p.dt:p.tpert+p.tau)/s.tg(end);
    [s.t_short,s.tp_short,s.x_short,s.xp_short]=thalamic_ode(p);
    
    
    s.t_short=s.t_short(p.tpert_idx:end);
    s.tp_short=s.tp_short(p.tpert_idx:end);
    s.x_short=s.x_short(p.tpert_idx:end,:);
    s.xp_short=s.xp_short(p.tpert_idx:end,:);
    
    
    p.perturbed=0;
    %% Sensitivity Function
    p.s0=[1 0 0 1];
    numeric_sens=0;
    [s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s,numeric_sens);
    s.t_sensitivity=s.t_short;
    s.sensitivity=s.sensitivity_direct;
    %% Approximations
    dtheta=p.dp*sum(s.prc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
    dpsi=p.dp*sum(s.irc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
    dpsi2=p.dp*sum(s.irc2(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
    
    %%
    final_vals_theta(jj,:)=dtheta;
    final_vals_psi(jj,:)=dpsi;   
    final_vals_psi2(jj,:)=dpsi2;
    set(er,'CData',final_vals_theta);
    pause(0.0001)
end



