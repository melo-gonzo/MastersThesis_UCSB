a_vals=0:0.01:1;
tau_vals=0:0.01:1;
final_vals=zeros(length(a_vals),length(tau_vals));
kk=1;
total=length(1:length(tau_vals))*length(2:length(tau_vals));
for jj=1:length(a_vals)
    for ii=2:length(tau_vals)
        kk=kk+1;
        if mod(kk,10)==0
            disp(round(100*kk/total,3))
        end
        p.a=p.a(1);
        p.dp=a_vals(jj);
        tau=tau_vals(ii);
        
        p.tau=tau;
        p.tpert=0;
        p.perturbed=1;
%         p.dp=p.a(1)*0.5;
        p.a=[p.a p.a+p.dp p.a];
        [s.t_short,s.tp_short,s.x_short,s.xp_short]=hopf_ode(p);
        p.perturbed=0;
        [~,s.tp,~,s.xp]=hopf_ode(p);
        %% Sensitivity Function
        p.s0=[1 0 0 1];
        [s.sensitivity_direct, s.greens, s.greens_vector]=sensitivity_function(p,s);
%         [s.t_sensitivity, s.sensitivity]=variation_of_constants(p,s);
        s.t_sensitivity=s.t_short;
        s.sensitivity=s.sensitivity_direct;
        %% Coordinate Recovery
        s.xy_cr=cr_function(p,s);
        s.xy_recovered = s.x_short + s.xy_cr;
        %%
        final_vals(jj,ii)=sqrt(s.xy_recovered(end,1)^2+s.xy_recovered(end,2)^2);
    end 
end 