function dtheta=delta_theta(t,x_theta,p,s)
s.prc(1,1)=0;
if length(p.phi_span)==2
    p.phi_span=p.phi_span_um;
end 

[t_p,x_p]=ode_interp_forward(t,[p.phi_span' s.x_short...
    s.prc(p.tpert_idx:p.tau_idx,2:end)]);
%x_p=[x y prc_x prc_y]
dtheta=x_p(3,1)*x_p(1,1)+x_p(4,1)*x_p(2,1);
