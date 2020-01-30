% Figuring out sensitivity function calculation
close all
clear all
clc
set_default_plot
dt=0.005;
tspan=0:dt:20*pi;


[tsol,xsol]=ode45(@(t,x) [x(1)-x(2)-(x(1)+x(2))*(x(1)^2+x(2)^2);...
    x(1)+x(2)+(x(1)-x(2))*(x(1)^2+x(2)^2)], tspan, [1;0]);
%% Sensitivity Direct
[~,sd]=ode45(@(t,sd) sd_f(t,sd,tsol,xsol), tspan, [0;0]);

plot(sd)

function sxy_dot=sd_f(t,sd,tsol,xsol)
[~,x_p]=ode_interp_forward(t,[tsol, xsol]);
% J=[cos(t) sin(t);
%     sin(t) cos(t)];
x=x_p(1);
y=x_p(2);
J(1,1)=1+2*x*(-x-y)-(x^2+y^2);
J(1,2)=-1+2*y*(-x-y)-(x^2+y^2);
J(2,1)=1+2*x*(x-y)+(x^2+y^2);
J(2,2)=1+2*y*(x-y)-(x^2+y^2);


sxy_dot(1,1)=J(1,1)*sd(1,1)+J(1,2)*sd(2,1)+x_p(1,1);
sxy_dot(2,1)=J(2,1)*sd(1,1)+J(2,2)*sd(2,1)+x_p(2,1);
end


%% Sensitiviyt Numeric 
for k=1:length(tspan)-1
    [tg,greens_2]=ode45(@(t,greens)...
        k_ode(t,greens,tsol,xsol,a), flip(tspan(k):dt:tspan(k+1)) ,1);
    greens_vector_2(k,1)=greens_2(end);
end 
greens_vector_2=(cumprod([greens_vector_2(1);greens_vector_2]));
    
%%%%%
%%%%%
[~,greens_vector]=ode45(@(t,greens) k_ode(t,greens,tsol,xsol,a), 0:dt:tsol(end) ,1);
greens_vector=flip(greens_vector);
%%%%%
%%%%%

% 
% % greens_vector=nan(length(tspan),1);
% % idx=1;
% % tot=length(tspan);
% % for t_solve=tspan
% %     if mod(idx,100)==0
% %         disp(round(100*idx/tot,3))
% %     end 
% %     try
% %     [~,greens]=ode45(@(t,greens) k_ode(t,greens,tsol,xsol,a), 0:dt:t_solve ,1);
% %     greens_vector(idx,:)=greens(end,:);
% %     catch
% %         greens_vector(idx,:)=1;
% %     end
% %     idx=idx+1;
% % end 
% % greens_vector=flip(greens_vector);
% 
% 
% [~,sn]=ode45(@(t,sdot)...
%     variation_of_constants(t,sdot,tsol,xsol,greens_vector), tspan, 0);
% %% 
% plot(tspan,sa);hold on;
% plot(tspan,sd);
% plot(tspan,sn/2);
% % plot(tspan,sn.*greens_vector)
% legend('$S_{analytic}$','$S_{direct}$','$S_{numerical}$')
% 
% % plot(tspan,greens_vector);
% %% Functions 
% function sd_dot=sd_f(t,sd,tsol,xsol,a,~)
% [~,x_p]=ode_interp_forward(t,[tsol,(xsol.*cos(tsol))]);
% dfdp=x_p;
% J=a*cos(t);
% sd_dot=J*sd+dfdp;
% end
% 
% function Gamma = k_ode(t,greens,~,~,a)
%     J=-a*cos(t);
%     Gamma(1,1)=greens(1)*J;
% end 
% 
% function sensitivity = variation_of_constants(t,~,tsol,xsol,greens_vector)
% %     [~,x_p]=ode_interp_forward(t,[tsol [greens_vector xsol.*cos(tsol)./greens_vector]]);
% %      sensitivity=x_p(2);
% 
%     [~,x_p]=ode_interp_forward(t,[tsol [greens_vector greens_vector.*xsol.*cos(tsol)./flip(greens_vector)]]);
%      sensitivity=x_p(2);
% 
% end 
% 
