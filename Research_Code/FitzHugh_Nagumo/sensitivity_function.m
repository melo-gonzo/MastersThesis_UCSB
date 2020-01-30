function [sensitivity_direct, greens, greens_vector]=sensitivity_function(p,s,numeric)
if numeric==1
disp('Calculating Greens Function')
end 

[~,sensitivity_direct]=ode45(@(t,sxy) sens_direct(t,sxy,p,s),...
    (p.tpert:p.dt:p.tpert+p.tau),[0,0]);
greens=nan;greens_vector=nan;

% options=odeset('RelTol',3e-10,'AbsTol',1e-10*ones(length(p.s0),1));
if numeric==1
    greens_vector=nan(length(p.tpert:p.dt:p.tpert+p.tau),4);
    idx=1;
    for t_solve=p.tpert:p.dt:(p.tpert+p.tau)
        try
            [t_greens,greens]=ode45(@(t,greens) gamma_ode(t,greens,p,s),...
                flip(p.tpert:p.dt:t_solve),p.s0);
            greens=flip(greens);t_greens=flip(t_greens);
            greens_vector(idx,:)=greens(1,:);
        catch
            greens_vector(idx,:)=[1,0,0,1];
        end
        idx=idx+1;
    end
    aaaa=1;
end 
    function [Gamma] = gamma_ode(t,greens,p,s)
    [~,x_p]=ode_interp_backward(t,[s.tp_short s.xp_short]);
    [J,~]=DF_nhopf(x_p,p);
    Gamma(1,1)=-(greens(1,1)*J(1,1)+greens(2,1)*J(2,1));
    Gamma(2,1)=-(greens(1,1)*J(1,2)+greens(2,1)*J(2,2));
    Gamma(3,1)=-(greens(3,1)*J(1,1)+greens(4,1)*J(2,1));
    Gamma(4,1)=-(greens(3,1)*J(1,2)+greens(4,1)*J(2,2));
    end
end

