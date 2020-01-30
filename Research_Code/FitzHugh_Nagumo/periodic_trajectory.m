function [tg,xg]=periodic_trajectory(t,x)
n=2;
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
