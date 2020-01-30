% Figuring out sensitivity function calculation
close all
clear all
clc
set_default_plot
a=1;
dt=0.001;
tspan=0:dt:4*pi;
[tsol,xsol]=ode45(@(t,x) a*x*cos(t), tspan, 1);
%% Analytic Sensitivity
sa=sin(tspan).*exp(a*sin(tspan));

%% Sensitivity Direct
[~,sd]=ode45(@(t,sd) sd_f(t,sd,tsol,xsol,a), tspan, 0);

%% Sensitiviyt Numeric 

greens_vector=nan(length(tspan),1);
idx=1;
% tspan_f=flip(tspan);
tot=length(tspan);
for t_solve=tspan(end)
%     if mod(idx,100)==0
%         disp(round(100*idx/tot,3))
%     end 
%     try
    [~,greens]=ode45(@(t,greens)...
        k_ode(t,greens,tsol,xsol,a), 0:dt:t_solve ,1);
%     greens_vector(idx,:)=greens(end,:);
%     catch
%         greens_vector(idx,:)=1;
%     end
%     idx=idx+1;
%     if idx==round(tot/2)
%         aaa=1;
%     end
end 
greens_vector=flip(greens);
% greens_vector=-cos(tsol);
% greens_vector=flip(greens_vector);

[~,sensitivity]=ode45(@(t,sdot)...
    variation_of_constants(t,sdot,tsol,xsol,greens_vector), tspan, 0);
%% 
plot(tspan,sa);hold on;
% plot(tspan,sd);
plot(tspan,sensitivity)
legend('S_a','S_n')

% plot(tspan,greens_vector);
%% Functions 
function sd_dot=sd_f(t,sd,tsol,xsol,a)
[~,x_p]=ode_interp_forward(t,[tsol,xsol.*cos(tsol)]);
dfdp=x_p;
J=a*cos(t);
sd_dot=J*sd+dfdp;
end

function Gamma = k_ode(t,greens,tsol,xsol,a)
    J=a*cos(t);
    Gamma(1,1)=-greens(1)*J;
end 

function sensitivity = variation_of_constants(t,sdot,tsol,xsol,greens_vector)
    [~,x_p]=ode_interp_forward(t,[tsol [greens_vector xsol.*cos(tsol)]]);
    sensitivity=x_p(1)*x_p(2);
end 

