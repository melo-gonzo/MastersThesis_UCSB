function [t_sensitivity,sensitivity]=variation_of_constants_2(p,s)
disp('Calculating Sensitivity Function')
% options=odeset('RelTol',3e-10,'AbsTol',1e-10*ones(length(p.x0),1));
[t_sensitivity,sensitivity]=ode45(@(t,S) VOC_ode(t,S,p,s),...
    p.tpert:p.dt:(p.tau+p.tpert), [0,0]);
    function Sdot=VOC_ode(t,S,p,s)
        xg=[pinv(s.greens_vector)'...
            s.x_short(length(0:p.dt:p.tpert):end,:)];
        
        
        [~,x_p]=ode_interp_forward(t,[(p.tpert:p.dt:p.tpert+p.tau)' xg]);
        
        
        Sdot(1,1)=x_p(1)*x_p(5)+x_p(2)*x_p(6);
        Sdot(2,1)=x_p(3)*x_p(5)+x_p(4)*x_p(6);
        
        
    end
end