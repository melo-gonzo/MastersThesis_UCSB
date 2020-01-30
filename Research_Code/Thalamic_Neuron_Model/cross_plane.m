function [cross,isterminal,direction] = cross_plane(~,f,prc_ic,f0)
    cross=prc_ic(1)*(f(1)-f0(1))+prc_ic(2)*(f(2)-f0(2))+prc_ic(3)*(f(3)-f0(3));
    isterminal=0;
    direction=1;
end