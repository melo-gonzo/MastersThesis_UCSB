function xy_cr = cr_function(p,s,sim)
if sim==0
% disp('Calculating Coordinate Recovery')
end 
phis=p.phi_span;
if sim==1
    p.phi_span_um=p.phi_span;
    p.phi_span=[p.phi_span(1) p.phi_span(end)];
    phis=phis(end);
end

p.dt_psi=mean(diff(s.prc(:,1)));

dtheta=p.dp*sum(s.prc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
dpsi=p.dp*sum(s.irc(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);
dpsi2=p.dp*sum(s.irc2(p.tpert_idx:p.tau_idx,2:end).*s.sensitivity,2);

ttheta=p.phi_span;
tpsi=p.phi_span;
xy_cr=nan(length(phis),3);


for k=1:length(phis)
	d_theta = dtheta(find(ttheta>=phis(k),1));
    d_psi = dpsi(find(tpsi>=phis(k),1));
    d_psi2 = dpsi2(find(tpsi>=phis(k),1));
    s1=find(s.prc(:,1)>=phis(k),1);
    s2=find(s.irc(:,1)>=phis(k),1);
    s3=find(s.irc2(:,1)>=phis(k),1);
    JAPR=[s.prc(s1,2) s.prc(s1,3) s.prc(s1,4);
        s.irc(s2,2) s.irc(s2,3) s.irc(s2,4);
        s.irc2(s3,2) s.irc2(s3,3) s.irc2(s3,4)];
    if k==length(phis)
        d_theta=dtheta(end);d_psi=dpsi(end);d_psi2=dpsi2(end);
    end
    xy_theory=transpose((JAPR)\[d_theta;d_psi;d_psi2]);
    xy_cr(k,:)=[xy_theory];
end 

% for k=1:length(phis)
% 	d_theta = dtheta(find(ttheta>=phis(k),1));
%     d_psi = dpsi(find(tpsi>=phis(k),1));
%     s1=find(s.prc(:,1)>=phis(k),1);
%     s2=find(s.irc(:,1)>=phis(k),1);
%     JAPR=[s.prc(s1,2) s.prc(s1,3) s.prc(s1,4);
%         s.irc(s2,2) s.irc(s2,3) s.irc(s2,4)];
%     if k==length(phis)
%         d_theta=dtheta(end);d_psi=dpsi(end);
%     end
%     xy_theory=transpose((JAPR)\[d_theta;d_psi]);
%     xy_cr(k,:)=[xy_theory];
% end 



