function [prc_val]=prc_ode(t,prc_sol,xg_f,p,s)
[t_p, x_p]=ode_interp_backward(t,[s.tg xg_f]);
[D,~]=DF_nhopf(x_p,p);
prc_val=-D'*prc_sol;
