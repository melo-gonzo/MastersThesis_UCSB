function sxy_dot=sens_direct(t,sxy,p,s)
[~,x_p]=ode_interp_forward(t,[s.tp_short s.x_short]);
[J,~]=DF_thalamic(x_p,p);
sxy_dot(1,1)=J(1,1)*sxy(1,1)+J(1,2)*sxy(2,1)+J(1,3)*sxy(3,1)+(1/p.cm);
sxy_dot(2,1)=J(2,1)*sxy(1,1)+J(2,2)*sxy(2,1)+J(2,3)*sxy(3,1)+0;
sxy_dot(3,1)=J(3,1)*sxy(1,1)+J(3,2)*sxy(2,1)+J(3,3)*sxy(3,1)+0;


% sxy_dot(1,1)=J(1,1)*sxy(1,1)+J(1,2)*sxy(2,1)+x_p(1,1);
% sxy_dot(2,1)=J(2,1)*sxy(1,1)+J(2,2)*sxy(2,1)+x_p(2,1);
