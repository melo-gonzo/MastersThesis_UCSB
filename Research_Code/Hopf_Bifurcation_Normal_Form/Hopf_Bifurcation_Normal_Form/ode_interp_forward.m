function [t_p,x_p]=ode_interp_forward(t,tgxg)
t_g=round(tgxg(:,1),5);
xg=tgxg(:,2:end);

%first find which time from tgamma is closest to the current time of the
%solver. even through we give a fixed time step, matlab solves using
%variable time step method, and outputs where we specified.
t1=find(t_g>t);
if t>=t_g(end)
    for k=1:min(size(xg))
        t_p=t_g(end);
        x_p(k,1)=xg(end,k);
    end 
else
    %we'll need to interpolate values of the periodic trajectory to be at
    %the exact value of time we are interested in. 
   t_p=t_g(t1(1)-1:t1(1));
    for k=1:min(size(xg))
        x_p_range(:,k)=xg(t1(1)-1:t1(1),k);
        x_p(k,1)=interp1(t_p,x_p_range(:,end),t);
    end 
end
%D is the jacobian of F(x_gamma) at a point x,y.
% [J,F]=DF_nhopf(x_p,p);
% D=D;
end 
