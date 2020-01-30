

% a_vals=0:0.5:25;
% tau_vals=round(linspace(0,round(s.tg(end),3),length(a_vals)),3);

a_vals=1;
tau_vals=round(linspace((3*pi/2)*s.tg(end)/(2*pi),round(s.tg(end),3),100),3);
tau_vals=round(linspace(5.2*s.tg(end)/(2*pi),5.3*s.tg(end)/(2*pi),100),3);
tau_vals=[7.02 7.021]
final_vals=zeros(length(a_vals),length(tau_vals));
final_coords=cat(3,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals,final_vals);
%[x y z xp yp zp xr yr zr time_remains_close]
% er=pcolor(a_vals,tau_vals,final_vals);colormap jet;
% set(er,'EdgeColor','none')


pause(0.01)
kk=1;
total=length(1:length(tau_vals))*length(2:length(tau_vals));
% plot3(s.xg(:,1),s.xg(:,2),s.xg(:,3));hold on;view([0 30])
for ii=1:length(a_vals)
    pausehere=1;
    for jj=length(tau_vals):-1:1
        kk=kk+1;
        if mod(kk,10)==0
            disp(round(100*kk/total,3))
        end
        p.ib=p.ib(1);
        p.dp=a_vals(ii);
        s.prc=s.prc_og;
        s.irc=s.irc_og;
        s.irc2=s.irc_og2;
        tau=round(s.tg(end)*0.2,3);
        p.tau=tau;
        p.tpert=tau_vals(jj);
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

        %% Coordinate Recovery
        s.xy_cr=cr_function(p,s,0);
        s.xy_recovered = s.x_short + s.xy_cr;
        %%
        pd=(s.xy_recovered-s.xp_short)./s.xp_short;
        k=find(sum(pd>0.05,2),1);
        if isempty(k);k=length(pd);end 
        td=s.t_short(k)-s.tg(p.tpert_idx);
        %%
        figure();title(num2str(tau_vals(jj)));hold on
        plot3(s.xg(:,1),s.xg(:,2),s.xg(:,3))
        plot3(s.xp_short(:,1),s.xp_short(:,2),s.xp_short(:,3),'b-.')
        plot3(s.xy_recovered(:,1),s.xy_recovered(:,2),s.xy_recovered(:,3),'r.')
        plot3(s.xy_recovered(k,1),s.xy_recovered(k,2),s.xy_recovered(k,3),'r*','markersize',10)
        view([0 30])
        pause(1)
%         del_last_line(2)
        %%
%         final_vals(jj,ii)=sqrt((s.xy_recovered(end,1)-s.xp_short(end,1))^2+(s.xy_recovered(end,2)-s.xp_short(end,2))^2+(s.xy_recovered(end,3)-s.xp_short(end,3))^2);
        final_vals(jj,ii)=td;
        final_coords(jj,ii,:)=[s.x_short(end,:) s.xp_short(end,:) s.xy_recovered(end,:) td];
%         set(er,'CData',final_vals);
        pause(0.0001)
    end
end
p.ib=p.ib(1);
% save thal_01 final_coords


