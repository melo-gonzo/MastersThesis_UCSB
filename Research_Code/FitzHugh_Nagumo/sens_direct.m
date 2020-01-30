function sxy_dot=sens_direct(t,sxy,p,s)
[~,x_p]=ode_interp_forward(t,[s.tp_short s.x_short]);
[J,~]=DF_nfhn(x_p,p);
sxy_dot(1,1)=J(1,1)*sxy(1,1)+J(1,2)*sxy(2,1)+0;
sxy_dot(2,1)=J(2,1)*sxy(1,1)+J(2,2)*sxy(2,1)+p.d(1);


% sxy_dot(1,1)=J(1,1)*sxy(1,1)+J(1,2)*sxy(2,1)+x_p(1,1);
% sxy_dot(2,1)=J(2,1)*sxy(1,1)+J(2,2)*sxy(2,1)+x_p(2,1);
