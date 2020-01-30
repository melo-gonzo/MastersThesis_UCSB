function [sensitivity_direct, greens, greens_vector, greens_rs]=sensitivity_function_2(p,s,numeric)
if numeric==1
    disp('Calculating Greens Function')
end

options=odeset('RelTol',3e-7,'AbsTol',1e-7*ones(2,1));
[~,sensitivity_direct]=ode45(@(t,sxy) sens_direct(t,sxy,p,s),...
    (p.tpert:p.dt:p.tpert+p.tau),[0,0],options);
greens=nan;greens_vector=nan;greens_rs=nan;

if numeric==1
    
    options=odeset('RelTol',3e-7,'AbsTol',1e-7*ones(4,1));
    [~,greens]=ode45(@(t,greens) gamma_ode(t,greens,p,s),...
        p.tpert:p.dt:p.tpert+p.tau ,p.s0, options);
    
    greens_rs=reshape(greens',2,2,length(p.tpert:p.dt:p.tpert+p.tau));
    
    gv=nan(2,2,length(length(p.tpert:p.dt:p.tpert+p.tau)));
    for k=1:length(p.tpert:p.dt:p.tpert+p.tau)
        %gv(:,:,k)=greens_rs(:,:,k)*inv(greens_rs(:,:,end-k+1));
%         gv(:,:,k)=inv(greens_rs(:,:,end-k+1));
        gv(:,:,k)=inv(greens_rs(:,:,k));
    end
    
    gv_rs=reshape(gv,4,length(p.tpert:p.dt:p.tpert+p.tau))';
    greens_vector=gv_rs;
end
% reshape(greens,4,length(p.tpert:p.dt:p.tpert+p.tau))' to get back

    function [Gamma] = gamma_ode(t,greens,p,s)
        [~,x_p]=ode_interp_backward(t,[s.tp_short s.xp_short]);
        [J,~]=DF_nhopf(x_p,p);
        Gamma(1,1)=(greens(1,1)*J(1,1)+greens(2,1)*J(2,1));
        Gamma(2,1)=(greens(1,1)*J(1,2)+greens(2,1)*J(2,2));
        Gamma(3,1)=(greens(3,1)*J(1,1)+greens(4,1)*J(2,1));
        Gamma(4,1)=(greens(3,1)*J(1,2)+greens(4,1)*J(2,2));
    end
end

