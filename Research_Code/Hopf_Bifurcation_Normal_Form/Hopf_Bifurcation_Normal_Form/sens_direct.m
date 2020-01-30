function sxy_dot=sens_direct(t,sxy,p,s)
%s'=J(t)*s+(df/dp)(t)
[~,x_p]=ode_interp_forward(t,[s.tp_short s.x_short]);
[J,~]=DF_nhopf(x_p,p);
sxy_dot(1,1)=J(1,1)*sxy(1,1)+J(1,2)*sxy(2,1)+x_p(1,1); 
sxy_dot(2,1)=J(2,1)*sxy(1,1)+J(2,2)*sxy(2,1)+x_p(2,1);


% t=0:0.0011:s.tp_short(end);
% DF=nan(length(t),4);
% for k=1:length(t)
%     [~,x_p]=ode_interp_forward(t(k),[s.tp_short s.x_short]);
%     [J,~]=DF_nhopf(x_p,p);
%     DF(k,:)=reshape(J,[1,4]);
% end 