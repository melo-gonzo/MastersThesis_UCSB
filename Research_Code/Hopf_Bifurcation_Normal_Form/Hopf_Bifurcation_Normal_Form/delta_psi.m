function dpsi=delta_psi(t,x,p,s)
if length(p.phi_span)==2
    p.phi_span=p.phi_span_um;
end 
[~,x_p]=ode_interp_forward(t,[p.phi_span' s.x_short...
    s.irc(p.tpert_idx:p.tau_idx,2:end)]);
% x_p=[x y irc_x irc_y];
[~,sens]=ode_interp_forward(t,[s.t_sensitivity*2*pi/s.tg(end),s.sensitivity]);
%sens = [sx sy];
dpsi=x_p(3,1)*(s.kappa(1)*sens(1)+x_p(1))+...
    x_p(4,1)*(s.kappa(1)*sens(2)+x_p(2));
