function xy_cr = cr_function(p,s,sim)
if sim==0
disp('Calculating Coordinate Recovery')
end 
T=s.tg(end);
% phis=mod(p.phi_span,T);
phis=p.phi_span;
if sim==1
    p.phi_span_um=p.phi_span;
    p.phi_span=[p.phi_span(1) p.phi_span(end)];
    phis=phis(end);
end

% [ttheta,dtheta]=ode45(@(t,x_theta)
% delta_theta(t,x_theta,p,s),p.phi_span,0);
% [tpsi,dpsi]=ode45(@(t,x_psi) delta_psi(t,x_psi,p,s),p.phi_span,0);
% dtheta=p.dp*dtheta;
% dpsi=p.dp*dpsi;
% p.dt_psi=(p.b(1)+p.d(1)*sqrt(-p.a(1)/p.c(1))^2)*p.dt;

p.dt_psi=mean(diff(s.prc(:,1)));

%%%
% dtheta=p.dp*cumsum(p.dt_psi*sum((s.prc(p.tpert_idx:p.tau_idx,2:end)).*[0 p.d(1)],2));
% dpsi=p.dp*cumsum(p.dt_psi*sum(s.irc(p.tpert_idx:p.tau_idx,2:end).*(s.kappa(1)*s.sensitivity+[0 p.d(1)]),2));
%%%
%%%
dtheta=p.dp*sum(s.prc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
dpsi=(p.dp*sum(s.irc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2));
%%%

ttheta=p.phi_span;
tpsi=p.phi_span;
xy_cr=nan(length(phis),2);
for k=1:length(phis)
	d_theta = dtheta(find(ttheta>=phis(k),1));
    d_psi = dpsi(find(tpsi>=phis(k),1));
    s1=find(s.prc(:,1)>=phis(k),1);
    s2=find(s.irc(:,1)>=phis(k),1);
    JAPR=[s.prc(s1,2) s.prc(s1,3);
        s.irc(s2,2) s.irc(s2,3)];
    if k==length(phis)
        d_theta=dtheta(end);d_psi=dpsi(end);
    end
    xy_theory=transpose((JAPR)\[d_theta;d_psi]);
    xy_cr(k,:)=[xy_theory];
end 
