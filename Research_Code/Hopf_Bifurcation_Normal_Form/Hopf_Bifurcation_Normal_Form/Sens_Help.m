% Figuring out sensitivity function calculation
close all
clear all
clc
set_default_plot
a=1;
dt=0.005;
tspan=0:dt:2*pi;
[tsol,xsol]=ode45(@(t,x) a*x*cos(t), tspan, 1);
%% Analytic Sensitivity
sa=sin(tspan).*exp(a*sin(tspan));

%% Sensitivity Direct
[~,sd]=ode45(@(t,sd) sd_f(t,sd,tsol,xsol,a), tspan, 0);
plot(sd)
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


% greens_vector=nan(length(tspan),1);
% idx=1;
% tot=length(tspan);
% for t_solve=tspan
%     if mod(idx,100)==0
%         disp(round(100*idx/tot,3))
%     end 
%     try
%     [~,greens]=ode45(@(t,greens) k_ode(t,greens,tsol,xsol,a), flip(0:dt:t_solve) ,1);
%     greens_vector(idx,:)=greens(end,:);
%     catch
%         greens_vector(idx,:)=1;
%     end
%     idx=idx+1;
% end 
% greens_vector=flip(greens_vector);


[~,sn]=ode45(@(t,sdot)...
    variation_of_constants(t,sdot,tsol,xsol,greens_vector), tspan, 0);
%% 
plot(tspan,sa);hold on;
plot(tspan,sd);
plot(tspan,sn/2);
% plot(tspan,sn.*greens_vector)
legend('$S_{analytic}$','$S_{direct}$','$S_{numerical}$')

% plot(tspan,greens_vector);
%% Functions 
function sd_dot=sd_f(t,sd,tsol,xsol,a,~)
[~,x_p]=ode_interp_forward(t,[tsol,(xsol.*cos(tsol))]);
dfdp=x_p;
J=a*cos(t);
sd_dot=J*sd+dfdp;
end

function Gamma = k_ode(t,greens,~,~,a)
    J=-a*cos(t);
    Gamma(1,1)=greens(1)*J;
end 

function sensitivity = variation_of_constants(t,~,tsol,xsol,greens_vector)
%     [~,x_p]=ode_interp_forward(t,[tsol [greens_vector xsol.*cos(tsol)./greens_vector]]);
%      sensitivity=x_p(2);

    [~,x_p]=ode_interp_forward(t,[tsol [greens_vector greens_vector.*xsol.*cos(tsol)./flip(greens_vector)]]);
     sensitivity=x_p(2);

end 

