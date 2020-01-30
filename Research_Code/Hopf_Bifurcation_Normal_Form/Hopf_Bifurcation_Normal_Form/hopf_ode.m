function [t,tp,xsol,xsolp]=hopf_ode(p)
tp=[];
tspan=p.tspan;
if p.perturbed~=0
    p.tspan=0:p.dt:round(p.tpert+p.tau,10);
    tspan=p.tspan;
end 
xsolp=[];
options=odeset('RelTol',3e-14,'AbsTol',1e-15*ones(length(p.x0),1));
if length(p.a)>1 || length(p.b)>1 || length(p.c)>1 || length(p.d)>1
    if isfield(p,'tpert')
        if ~isfield(p,'tau')
            p.tau=1;
        end 
        p.tpert=[0 p.tpert(1) p.tpert+p.tau p.tspan(end)];
    else
        p.tau=1;
        p.tpert=0:p.tau:max([length(p.a),length(p.b),length(p.c),length(p.d)])-1;
        p.tpert=[0 p.tpert];
        p.tpert(end)=p.tspan(end);
    end
    segs=ones(4,length(p.tpert)-1);
    segs(1,:)=p.a.*segs(1,:);
    segs(2,:)=p.b.*segs(2,:);
    segs(3,:)=p.c.*segs(3,:);
    segs(4,:)=p.d.*segs(4,:);
else
    segs=[p.a;p.b;p.c;p.d];
    p.tpert=nan;
end



p.tpert=round(p.tpert,3);
p.a=segs(1,1);
p.b=segs(2,1);
p.c=segs(3,1);
p.d=segs(4,1);
if length(p.tspan)>1
[t,xsol]=ode45(@(t,x) hopf_ode_equations(t,x,p),p.tspan,p.x0,options);
else
    [t,xsol]=deal(nan);
end 
tidxhit1=[];
tidxhit2=[];
if length(p.tpert)<1
    [tseg,xseg]=deal(nan);
else
    for k=2:length(p.tpert)
        p.a=segs(1,k-1);
        p.b=segs(2,k-1);
        p.c=segs(3,k-1);
        p.d=segs(4,k-1);
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
            [tseg,xseg]=ode45(@(t,x) hopf_ode_equations(t,x,p),p.tspan,p.x0,options);
            tp=[tp;tseg];
            xsolp=[xsolp;xseg];
        elseif length(p.tspan)==0
            [tseg,xseg]=deal(nan);
        else
        end
    end
end 








