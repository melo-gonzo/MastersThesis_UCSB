function [xg,tg,locs]=xg_function(p,numeric)
% if numeric==1
disp('Calculating xg/tg')
% options=odeset('RelTol',3e-14,'AbsTol',1e-18*ones(length(p.x0),1));
[t,~,x,~]=thalamic_ode(p);
n=3;
[~,locs]=findpeaks(x(:,1));
if length(locs)<1
    xg=x;
    tg=t;
    disp('Run Longer Simulation')
elseif length(locs)==1
    tg=t(1:locs(end))-t(1);
    xg=nan(length(tg),n);
    for ii=1:n
        xg(:,ii)=x(1:locs(end),ii);
    end
else
    tg=t(locs(end-1):locs(end))-t(locs(end-1));
    xg=nan(length(tg),n);
    for ii=1:n
        xg(:,ii)=x(locs(end-1):locs(end),ii);
    end
end
locs=[locs(end-1) locs(end)];
% else 
%     rr=sqrt(-p.a(1)/p.c(1));
%     T=round(2*pi/(p.b(1)+p.d(1)*(rr^2)),3);
%     t=0:2*pi/((T/p.dt)-1):2*pi;
%     xg=[(rr*cos(t))' (rr*sin(t))'];
%     tg=0:p.dt:T;
% end 