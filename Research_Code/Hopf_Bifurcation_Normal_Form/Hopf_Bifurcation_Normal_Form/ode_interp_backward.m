function [t_p, x_p]=ode_interp_backward(t,tgxg)
% returns the time and function points at the specif time of the ode
% solver
% t_g=flip(tg);
tgxg=flip(tgxg);
t_g=tgxg(:,1);
xg=tgxg(:,2:end);
%first find which time from tgamma is closest to the current time of the
%solver. even through we give a fixed time step, matlab solves using
%variable time step method, and outputs where we specified.
t1=find(t_g<t);
if t==0
    for k=1:min(size(xg))
        t_p=0;
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
% [D,F]=DF_nhopf(x_p,p);
% D=D';
end 
