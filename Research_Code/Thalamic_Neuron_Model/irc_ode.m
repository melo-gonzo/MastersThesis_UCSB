function [irc_val]=irc_ode(t,irc_sol,xg_f,p,s)
[~,x_p]=ode_interp_backward(t,[s.tg xg_f]);
[D,~]=DF_thalamic(x_p,p);
irc_val=(s.kappa(s.scale_vector_idx)*eye(3)-D')*irc_sol;
