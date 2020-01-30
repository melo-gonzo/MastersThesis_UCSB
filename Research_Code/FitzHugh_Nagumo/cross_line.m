function [value,isterminal,direction]=cross_line(~,x,prc_ic,f0)
%     global f0 
%     value = x(2);
    value=prc_ic(1)*(x(1)-f0(1))+prc_ic(2)*(x(2)-f0(2));
    isterminal=0;
    direction=0;
end