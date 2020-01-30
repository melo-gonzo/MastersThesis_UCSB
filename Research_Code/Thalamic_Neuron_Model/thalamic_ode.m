function [t,tp,xsol,xsolp]=thalamic_ode(p)
tp=[];
tspan=p.tspan;
if p.perturbed~=0
    p.tspan=0:p.dt:round(p.tpert+p.tau,10);
    tspan=p.tspan;
end 
xsolp=[];
options=odeset('RelTol',3e-14,'AbsTol',1e-18*ones(length(p.x0),1));
if length(p.ib)>1
    if isfield(p,'tpert')
        if ~isfield(p,'tau')
            p.tau=1;
        end 
        p.tpert=[0 p.tpert(1) p.tpert+p.tau p.tspan(end)];
    else
        p.tau=1;
        p.tpert=0:p.tau:max([length(p.ib)])-1;
        p.tpert=[0 p.tpert];
        p.tpert(end)=p.tspan(end);
    end
    segs=ones(1,length(p.tpert)-1);
    segs(1,:)=p.ib.*segs(1,:);
else
    segs=[p.ib];
    p.tpert=nan;
end



p.tpert=round(p.tpert,3);
p.ib=segs(1,1);
if length(p.tspan)>1
[t,xsol]=ode45(@(t,x) thalamic_ode_equations(t,x,p),p.tspan,p.x0,options);
else
    [t,xsol]=deal(nan);
end 
tidxhit1=[];
tidxhit2=[];
if length(p.tpert)<1
    [tseg,xseg]=deal(nan);
else
    for k=2:length(p.tpert)
        p.ib=segs(1,k-1);
        %%%
        tidx1=find(p.tpert(k-1)<=tspan,1);
        tidx2=find(p.tpert(k)<=tspan,1);
        if ismember(tidx1,tidxhit2) & tidx1~=1
            tidx1=tidx1+1;
        end 
        tidxhit1(end+1)=tidx1;
        tidxhit2(end+1)=tidx2;
        p.tspan=nan;
        if tidx1<=length(tspan) | tidx2<=length(tspan)
            p.tspan=tspan(tidx1:tidx2);
        end
        if isempty(xsolp) | isnan(xsolp)
            p.x0=p.x0;
        else
            p.x0=xsolp(end,:);
        end
        if length(p.tspan)>1
            [tseg,xseg]=ode45(@(t,x) thalamic_ode_equations(t,x,p),p.tspan,p.x0,options);
            tp=[tp;tseg];
            xsolp=[xsolp;xseg];
        elseif length(p.tspan)==0
            [tseg,xseg]=deal(nan);
        else
        end
    end
end 








